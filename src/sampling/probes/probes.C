/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
08-06-2020  Jeff Heylmun  : Reload exisiting probe files and move probes to
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

#include "probes.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "IOmanip.H"
#include "mapPolyMesh.H"
#include "polyPatch.H"
#include "SortableList.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(probes, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        probes,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::probes::findElements(const fvMesh& mesh, const bool print)
{
    if (debug)
    {
        Info<< "probes: resetting sample locations" << endl;
    }

    elementList_.clear();
    elementList_.setSize(size());

    faceList_.clear();
    faceList_.setSize(size());

    forAll(*this, probei)
    {
        const vector& location = origPoints_[probei];

        const label celli = mesh.findCell(location);

        elementList_[probei] = celli;
        faceList_[probei] = findFaceIndex(mesh, celli, location);

        if (debug && (elementList_[probei] != -1 || faceList_[probei] != -1))
        {
            Pout<< "probes : found point " << location
                << " in cell " << elementList_[probei]
                << " and face " << faceList_[probei] << endl;
        }
    }


    boolList foundList(size(), true);
    label nBadProbes = 0;

    // Check if all probes have been found.
    forAll(elementList_, probei)
    {
        const vector& location = origPoints_[probei];

        if (Pstream::parRun())
        {
            labelList celliL(Pstream::nProcs(), elementList_[probei]);
            labelList faceiL(Pstream::nProcs(), faceList_[probei]);

            // Check at least one processor with cell.
            Pstream::gatherList(celliL);
            Pstream::scatterList(celliL);
            Pstream::gatherList(faceiL);
            Pstream::scatterList(faceiL);

            //- Find first actual cell and face
            bool set = false;
            forAll(celliL, proci)
            {
                if
                (
                    (celliL[proci] >= 0 && faceiL[proci] >= 0)
                 && celliL[proci] == elementList_[probei]
                 && faceiL[proci] == faceList_[probei]
                 && !set
                )
                {
                    set = true;
                }
                else
                {
                    elementList_[probei] = -1;
                    faceList_[probei] = -1;
                }

            }
        }
        label celli = elementList_[probei];
        label facei = faceList_[probei];

        reduce(celli, maxOp<label>());
        reduce(facei, maxOp<label>());

        if (celli == -1)
        {
            foundList[probei] = false;

            if (Pstream::master() && print)
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
        reduce(foundList[probei], andOp<bool>());
        if (!foundList[probei])
        {
            nBadProbes++;
        }
    }

    if (!returnReduce(mesh.nFaces(), sumOp<label>()))
    {
        return;
    }

    if (adjustLocations_ && nBadProbes > 0)
    {
        if (print)
        {
            Info<< nl
                << nBadProbes << " probes were not found in any domain." << nl
                << "These probes are being moved to the nearest patch face." << nl
                << endl;
        }

        forAll(foundList, probei)
        {
            if (foundList[probei])
            {
                continue;
            }

            scalarField distance(mag(origPoints_[probei] - mesh.cellCentres()));
            label celli(findMin(distance));
            scalar minDistance = great;
            if (celli >= 0)
            {
                minDistance = distance[celli];
            }

            scalarList trueMinDistancel(Pstream::nProcs(), minDistance);
            Pstream::gatherList(trueMinDistancel);
            Pstream::scatterList(trueMinDistancel);
            label minI = max(findMin(trueMinDistancel), 0);

            if (Pstream::myProcNo() == minI)
            {
                elementList_[probei] = celli;
                faceList_[probei] = findFaceIndex(mesh, celli, origPoints_[probei]);
                operator[](probei) = mesh.faceCentres()[faceList_[probei]];

                if (print)
                {
                    Pout<< "Moved probe " << probei << nl
                        << "    Original position: " << origPoints_[probei] << nl
                        << "    New position: " << operator[](probei) << nl
                        << "    Located in cell " << elementList_[probei]
                        << ", face " << faceList_[probei] << endl;
                }
            }
            else
            {
                elementList_[probei] = -1;
                faceList_[probei] = -1;
            }
        }
        if (print)
        {
            Info<<endl;
        }
    }
}


Foam::label Foam::probes::findFaceIndex
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
            vector dist = mesh.faceCentres()[facei] - pt;
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


Foam::label Foam::probes::prepare()
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

                        // Do not overwrite files if the number of probes has
                        // changed
                        if (nOldProbes != size())
                        {
                            fileName oldProbeDir(probeDir);
                            probeDir = probeDir/".."/mesh_.time().timeName();
                            probeDir.clean();

                            if (Pstream::master())
                            {
                                WarningInFunction
                                    << "The number of probes in " << oldProbeDir
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
                fout<< "# Probe " << probei << ' ' << origPoints_[probei]
                    << endl;
            }

            fout<< '#' << setw(IOstream::defaultPrecision() + 6)
                << "Probe";

            forAll(*this, probei)
            {
                fout<< ' ' << setw(w) << probei;
            }
            fout<< endl;

            fout<< '#' << setw(IOstream::defaultPrecision() + 6)
                << "Time" << endl;

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

Foam::probes::probes
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
    fixedLocations_(true),
    interpolationScheme_("cell"),
    adjustLocations_(false),
    append_(false)
{
    read(dict);
}


Foam::probes::probes
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
    fixedLocations_(true),
    interpolationScheme_("cell"),
    adjustLocations_(false),
    append_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::probes::~probes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::probes::read(const dictionary& dict)
{
    dict.lookup("probeLocations") >> *this;
    origPoints_ = refCast<pointField>(*this);
    dict.lookup("fields") >> fieldSelection_;


    dict.readIfPresent("fixedLocations", fixedLocations_);
    if (dict.readIfPresent("interpolationScheme", interpolationScheme_))
    {
        if (!fixedLocations_ && interpolationScheme_ != "cell")
        {
            WarningInFunction
                << "Only cell interpolation can be applied when "
                << "not using fixedLocations.  InterpolationScheme "
                << "entry will be ignored";
        }
    }

    dict.readIfPresent("adjustLocations", adjustLocations_);
    dict.readIfPresent("append", append_);

    // Initialise cells to sample from supplied locations
    findElements(mesh_, true);

    prepare();

    return true;
}


bool Foam::probes::execute()
{
    return true;
}


bool Foam::probes::write()
{
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


void Foam::probes::updateMesh(const mapPolyMesh& mpm)
{
    DebugInfo<< "probes: updateMesh" << endl;

    if (&mpm.mesh() != &mesh_)
    {
        return;
    }

    if (fixedLocations_)
    {
        findElements(mesh_, false);
    }
    else
    {
        if (debug)
        {
            Info<< "probes: remapping sample locations" << endl;
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


void Foam::probes::movePoints(const polyMesh& mesh)
{
    DebugInfo<< "probes: movePoints" << endl;

    if (fixedLocations_ && &mesh == &mesh_)
    {
        findElements(mesh_, false);
    }
}


// ************************************************************************* //
