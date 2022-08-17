/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
08-06-2020  Synthetik Applied Technologies: |   Reload exisiting probe files and
                                                move blastProbes to
                                                the nearest cell/face.
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "blastProbes.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "IOmanip.H"
#include "mapPolyMesh.H"
#include "polyPatch.H"
#include "SortableList.H"
#include "IFstream.H"
#include "vtkWriteOps.H"
#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blastProbes, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        blastProbes,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::blastProbes::findElements
(
    const fvMesh& mesh,
    const bool print,
    const bool movePts
)
{
    if (debug)
    {
        Info<< "blastProbes: resetting sample locations" << endl;
    }

    elementList_.clear();
    elementList_.setSize(size());

    faceList_.clear();
    faceList_.setSize(size());

    boolList foundList(size(), false);
    label nBadProbes = 0;

    forAll(*this, probei)
    {
        const vector& location = operator[](probei);

        const label celli = mesh.findCell(location);

        elementList_[probei] = celli;
        faceList_[probei] = findFaceIndex(mesh, celli, location);
    }

    // Check if all blastProbes have been found.
    forAll(elementList_, probei)
    {
        const vector& location = operator[](probei);
        label celli = elementList_[probei];
        label facei = faceList_[probei];

        if (Pstream::parRun())
        {
            label proc = -1;
            // Favor keeping cells at the old cell centre if it is still valid
            if (celli >= 0 && facei >= 0)
            {
                if
                (
                    mag
                    (
                        mesh.cellCentres()[celli]
                      - elementLocations_[probei]
                    ) < small
                )
                {
                    proc = Pstream::myProcNo();
                }
            }
            if (returnReduce(proc, maxOp<label>()) < 0)
            {
                proc = (celli < 0 || facei < 0 ? -1 : Pstream::myProcNo());
            }
            reduce(proc, maxOp<label>());
            if (proc < 0 || proc != Pstream::myProcNo())
            {
                celli = -1;
                facei = -1;
                elementList_[probei] = -1;
                faceList_[probei] = -1;
                elementLocations_[probei] = vector(-great, -great, -great);
            }
            else
            {
                elementLocations_[probei] = mesh.cellCentres()[celli];
            }
            reduce(facei, maxOp<label>());
            reduce(celli, maxOp<label>());
            reduce(elementLocations_[probei], maxOp<vector>());
        }

        if ((elementList_[probei] != -1 || faceList_[probei] != -1))
        {
            foundList[probei] = true;
            if (debug)
            {
                Pout<< "blastProbes : found point " << location
                    << " in cell " << elementList_[probei]
                    << " cell centre " << mesh.cellCentres()[celli]
                    << " and face " << faceList_[probei] << endl;
            }
        }
        else if (celli == -1)
        {
            if (print)
            {
                WarningInFunction
                    << "Did not find location " << location
                    << " in any cell. Skipping location." << endl;
            }
        }
        else if (facei == -1)
        {
            if (print)
            {
                WarningInFunction
                    << "Did not find location " << location
                    << " in any face. Skipping location." << endl;
            }
        }
        else
        {
            // Make sure location not on two domains.
            if (elementList_[probei] != -1 && elementList_[probei] != celli)
            {
                if (print)
                {
                    WarningInFunction
                        << "Location " << location
                        << " seems to be on multiple domains:"
                        << " cell " << elementList_[probei]
                        << " on my domain " << Pstream::myProcNo()
                            << " and cell " << celli << " on some other domain."
                        << endl
                        << "This might happen if the probe location is on"
                        << " a processor patch. Change the location slightly"
                        << " to prevent this." << endl;
                }
            }

            if (faceList_[probei] != -1 && faceList_[probei] != facei)
            {
                if (print)
                {
                    WarningInFunction
                        << "Location " << location
                        << " seems to be on multiple domains:"
                        << " cell " << faceList_[probei]
                        << " on my domain " << Pstream::myProcNo()
                            << " and face " << facei << " on some other domain."
                        << endl
                        << "This might happen if the probe location is on"
                        << " a processor patch. Change the location slightly"
                        << " to prevent this." << endl;
                }
            }
        }
        reduce(foundList[probei], orOp<bool>());
        if (!foundList[probei])
        {
            nBadProbes++;
        }
    }

    needUpdate_ = false;

    if (!returnReduce(mesh.nFaces(), sumOp<label>()))
    {
        return;
    }

    if (nBadProbes > 0 && movePts)
    {
        if (print)
        {
            Info<< nl
                << nBadProbes << " blastProbes were not found in any domain." << nl
                << "These blastProbes are being moved to the nearest patch face."
                << nl
                << endl;
        }

        forAll(foundList, probei)
        {
            vector origPoint = operator[](probei);
            if (foundList[probei])
            {
                continue;
            }

            scalar minDistance = great;
            label facei = -1;
            for
            (
                label faceI = mesh.nInternalFaces();
                faceI < mesh.nFaces();
                faceI++
            )
            {
                scalar newDistance
                (
                    mag
                    (
                        origPoint
                      - mesh.faceCentres()[faceI]
                    )
                );
                if (newDistance < minDistance)
                {
                    facei = faceI;
                    minDistance = newDistance;
                }
            }

            scalar trueMinDistance =
                returnReduce(minDistance, minOp<scalar>());
            vector pt(-great, -great, -great);

            if (mag(trueMinDistance - minDistance) < small)
            {
                label patchi = mesh.boundaryMesh().whichPatch(facei);
                faceList_[probei] = facei;
                label localFaceID =
                    facei - mesh.boundaryMesh()[patchi].start();
                elementList_[probei] =
                    mesh.boundaryMesh()[patchi].faceCells()[localFaceID];
                pt = mesh.cellCentres()[elementList_[probei]];

                if (print || debug)
                {
                    Pout<< "Moved probe " << probei << nl
                        << "    Original position: " << origPoint << nl
                        << "    New position: " << pt << nl
                        << "    Located in cell " << elementList_[probei]
                        << ", face " << faceList_[probei] << endl;
                }
            }
            else
            {
                elementList_[probei] = -1;
                faceList_[probei] = -1;
            }
            reduce(pt, maxOp<vector>());
            operator[](probei) = pt;
            elementLocations_[probei] = pt;
        }
        if (print)
        {
            Info<<endl;
        }
    }
}


Foam::label Foam::blastProbes::findFaceIndex
(
    const fvMesh& mesh,
    const label celli,
    const vector& pt
) const
{
    if (celli != -1)
    {
        const labelList& cellFaces = mesh.cells()[celli];
        scalar minDistance = great;
        label minFaceID = -1;
        forAll(cellFaces, i)
        {
            label facei = cellFaces[i];
            const vector& cellCentre = mesh.cellCentres()[celli];
            vector dist = mesh.faceCentres()[facei] - cellCentre;
            if (mag(dist) < minDistance)
            {
                minDistance = mag(dist);
                minFaceID = facei;
            }
        }
        return minFaceID;
    }
    return -1;
}


Foam::label Foam::blastProbes::prepare()
{
    const label nFields = classifyFields();

    // adjust file streams
    if (Pstream::master())
    {
        wordHashSet currentFields;

        currentFields.insert(scalarFields_);
        currentFields.insert(vectorFields_);
        currentFields.insert(sphericalTensorFields_);
        currentFields.insert(symmTensorFields_);
        currentFields.insert(tensorFields_);

        currentFields.insert(surfaceScalarFields_);
        currentFields.insert(surfaceVectorFields_);
        currentFields.insert(surfaceSphericalTensorFields_);
        currentFields.insert(surfaceSymmTensorFields_);
        currentFields.insert(surfaceTensorFields_);

        if (debug)
        {
            Info<< "Probing fields: " << currentFields << nl
                << "Probing locations: " << *this << nl
                << endl;
        }


        fileName probeDir;
        fileName probeSubDir = name();

        if (mesh_.name() != polyMesh::defaultRegion)
        {
            probeSubDir = probeSubDir/mesh_.name();
        }
        probeSubDir = "postProcessing"/probeSubDir;

        if (Pstream::parRun())
        {
            // Put in undecomposed case
            // (Note: gives problems for distributed data running)
            probeDir = mesh_.time().path()/".."/probeSubDir;
        }
        else
        {
            probeDir = mesh_.time().path()/probeSubDir;
        }

        wordList times(readDir(probeDir, fileType::directory));

        // Sort times
        {
            SortableList<scalar> sTimes(times.size());
            forAll(sTimes, ti)
            {
                IStringStream is(times[ti]);
                sTimes[ti] = readScalar(is);
            }
            sTimes.sort();
            wordList oldTimes(times);
            forAll(sTimes, ti)
            {
                times[ti] = oldTimes[sTimes.indices()[ti]];
            }
        }

        word timeName;
        if (append_ && times.size())
        {
            timeName = times[0];
        }
        else
        {
            timeName = mesh_.time().timeName();
        }
        probeSubDir = probeSubDir/timeName;

        if (Pstream::parRun())
        {
            // Put in undecomposed case
            // (Note: gives problems for distributed data running)
            probeDir = mesh_.time().path()/".."/probeSubDir;
        }
        else
        {
            probeDir = mesh_.time().path()/probeSubDir;
        }
        // Remove ".."
        probeDir.clean();

        // ignore known fields, close streams for fields that no longer exist
        forAllIter(HashPtrTable<OFstream>, probeFilePtrs_, iter)
        {
            if (!currentFields.erase(iter.key()))
            {
                if (debug)
                {
                    Info<< "close probe stream: " << iter()->name() << endl;
                }

                delete probeFilePtrs_.remove(iter);
            }
        }

        // currentFields now just has the new fields - open streams for them
        forAllConstIter(wordHashSet, currentFields, iter)
        {
            const word& fieldName = iter.key();

            // Create directory if does not exist.
            mkDir(probeDir);

            // Read old file and store stream as a list of strings
            wordList oldValues;
            if
            (
                exists(fileName(probeDir/fieldName))
             && times[0] != mesh_.time().timeName()
             && append_
            )
            {
                label nOldProbes = 0;
                bool header = true;
                IFstream is(fileName(probeDir/fieldName));
                string line;

                while (is.good())
                {
                    is.getLine(line);

                    if (line[0] == '#')
                    {
                        nOldProbes++;
                    }
                    else if (header)
                    {
                        header = false;
                        nOldProbes -= 2;

                        // Do not overwrite files if the number of blastProbes has
                        // changed
                        if (nOldProbes != size())
                        {
                            fileName oldProbeDir(probeDir);
                            probeDir = probeDir/".."/mesh_.time().timeName();
                            probeDir.clean();

                            if (Pstream::master())
                            {
                                WarningInFunction
                                    << "The number of blastProbes in " << oldProbeDir
                                    << nl
                                    << "    is not the same as the previous file."
                                    << nl
                                    << "    The previous probe file will not be"
                                    << " overwritten. " << nl
                                    << "    Writing to "
                                    << probeDir << endl;
                            }
                            break;
                        }
                    }

                    if (!header)
                    {
                        oldValues.append(line);
                    }
                }
            }

            OFstream* fPtr = new OFstream(probeDir/fieldName);

            OFstream& fout = *fPtr;

            if (debug)
            {
                Info<< "open probe stream: " << fout.name() << endl;
            }

            probeFilePtrs_.insert(fieldName, fPtr);

            unsigned int w = IOstream::defaultPrecision() + 7;

            forAll(*this, probei)
            {
                fout<< "# Probe " << probei << ' ' << operator[](probei)
                    << endl;
            }


            fout<< '#' << setw(w) << "Time";
            forAll(*this, probei)
            {
                fout<< ' ' << setw(w) << probei;
            }
            fout<< endl;

            // Add old values to new output
            if (oldValues.size())
            {
                forAll(oldValues, i)
                {
                    IStringStream isLine(oldValues[i]);
                    scalar t = readScalar(isLine);

                    if (t <= mesh_.time().value())
                    {
                        fout << word(oldValues[i]) << nl;
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }
    }

    return nFields;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blastProbes::blastProbes
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    pointField(0),
    mesh_
    (
        refCast<const fvMesh>
        (
            t.lookupObject<objectRegistry>
            (
                dict.lookupOrDefault("region", polyMesh::defaultRegion)
            )
        )
    ),
    loadFromFiles_(false),
    fieldSelection_(),
    fixedLocations_(false),
    interpolationScheme_("cell"),
    append_(false)
{
    read(dict);
}


Foam::blastProbes::blastProbes
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObject(name),
    pointField(0),
    mesh_(refCast<const fvMesh>(obr)),
    loadFromFiles_(loadFromFiles),
    fieldSelection_(),
    fixedLocations_(false),
    interpolationScheme_("cell"),
    append_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blastProbes::~blastProbes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::blastProbes::read(const dictionary& dict)
{
    dict.lookup("probeLocations") >> *this;
    dict.lookup("fields") >> fieldSelection_;


    dict.readIfPresent("fixedLocations", fixedLocations_);
    if
    (
        dict.readIfPresent
        (
            "interpolationScheme",
            interpolationScheme_
        )
    )
    {
        if (!fixedLocations_ && interpolationScheme_ != "cell")
        {
            WarningInFunction
                << "Only cell interpolation can be applied when "
                << "not using fixedLocations.  InterpolationScheme "
                << "entry will be ignored";
        }
    }

    dict.readIfPresent("append", append_);
    if (!elementLocations_.size() || !fixedLocations_)
    {
        elementLocations_.clear();
        elementLocations_.setSize(size());
        elementLocations_ = Zero;

        // Initialise cells to sample from supplied locations
        findElements
        (
            mesh_,
            true,
            dict.lookupOrDefault("adjustLocations", false)
        );
    }
    prepare();

    Switch writeVTK(dict.lookupOrDefault("writeVTK", false));

    if (writeVTK && Pstream::master())
    {
        IOstream::streamFormat writeFormat = IOstream::ASCII;
        if (dict.found("writeFormat"))
        {
            writeFormat = IOstream::formatEnum
            (
                dict.lookup("writeFormat")
            );
        }
        const Time& time(mesh_.time());
        const fileName path = time.rootPath()/time.globalCaseName()/"VTK";
        mkDir(path);

        const fileName filePath = path/this->name() + ".vtk";

        bool binary(writeFormat == IOstream::BINARY);
        ofstream os(filePath, std::ios::binary);

        vtkWriteOps::writeHeader(os, binary, "sampleSurface");
        os << "DATASET POLYDATA" << nl;

        // Write vertex coords
        os  << "POINTS " << elementLocations_.size() << " float" << nl;

        List<floatScalar> po(elementLocations_.size()*3);
        label ind = 0;
        forAll(elementLocations_, pointi)
        {
            const point& pt = elementLocations_[pointi];
            forAll(pt, cmpt)
            {
                po[ind++] = float(pt[cmpt]);
            }
        }
        vtkWriteOps::write(os, binary, po);
    }

    return true;
}


bool Foam::blastProbes::execute()
{
    return true;
}


bool Foam::blastProbes::write()
{
    if (needUpdate_)
    {
        findElements(mesh_, true);
    }
    if (size() && prepare())
    {
        sampleAndWrite(scalarFields_);
        sampleAndWrite(vectorFields_);
        sampleAndWrite(sphericalTensorFields_);
        sampleAndWrite(symmTensorFields_);
        sampleAndWrite(tensorFields_);

        sampleAndWriteSurfaceFields(surfaceScalarFields_);
        sampleAndWriteSurfaceFields(surfaceVectorFields_);
        sampleAndWriteSurfaceFields(surfaceSphericalTensorFields_);
        sampleAndWriteSurfaceFields(surfaceSymmTensorFields_);
        sampleAndWriteSurfaceFields(surfaceTensorFields_);
    }

    return true;
}


void Foam::blastProbes::updateMesh(const mapPolyMesh& mpm)
{
    DebugInfo<< "blastProbes: updateMesh" << endl;

    if (&mpm.mesh() != &mesh_)
    {
        return;
    }

    if (!fixedLocations_)
    {
        needUpdate_ = true;
    }
    else
    {
        if (debug)
        {
            Info<< "blastProbes: remapping sample locations" << endl;
        }

        // 1. Update cells
        {
            DynamicList<label> elems(elementList_.size());

            const labelList& reverseMap = mpm.reverseCellMap();
            forAll(elementList_, i)
            {
                label celli = elementList_[i];
                label newCelli = reverseMap[celli];
                if (newCelli == -1)
                {
                    // cell removed
                }
                else if (newCelli < -1)
                {
                    // cell merged
                    elems.append(-newCelli - 2);
                }
                else
                {
                    // valid new cell
                    elems.append(newCelli);
                }
            }

            elementList_.transfer(elems);
        }

        // 2. Update faces
        {
            DynamicList<label> elems(faceList_.size());

            const labelList& reverseMap = mpm.reverseFaceMap();
            forAll(faceList_, i)
            {
                label facei = faceList_[i];
                label newFacei = reverseMap[facei];
                if (newFacei == -1)
                {
                    // face removed
                }
                else if (newFacei < -1)
                {
                    // face merged
                    elems.append(-newFacei - 2);
                }
                else
                {
                    // valid new face
                    elems.append(newFacei);
                }
            }

            faceList_.transfer(elems);
        }
    }
}


void Foam::blastProbes::movePoints(const polyMesh& mesh)
{
    DebugInfo<< "blastProbes: movePoints" << endl;

    if (!fixedLocations_ && &mesh == &mesh_)
    {
        needUpdate_ = true;
    }
}


// ************************************************************************* //
