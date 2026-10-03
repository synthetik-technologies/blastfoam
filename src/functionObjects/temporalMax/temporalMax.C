/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           |
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "temporalMax.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "Time.H"
#include "OFstream.H"
#include "Pstream.H"
#include "polyMeshMap.H"
#include "OSspecific.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(temporalMax, 0);
    addToRunTimeSelectionTable(functionObject, temporalMax, dictionary);

    // Aliases
    class temporalFieldMax : public temporalMax
    {
    public:
        TypeName("temporalFieldMax");
        temporalFieldMax(const word& name, const Time& runTime, const dictionary& dict)
        : temporalMax(name, runTime, dict) {}
    };
    defineTypeNameAndDebug(temporalFieldMax, 0);
    addToRunTimeSelectionTable(functionObject, temporalFieldMax, dictionary);

    class coarseMeshMaxAlias : public temporalMax
    {
    public:
        TypeName("coarseMeshMax");
        coarseMeshMaxAlias(const word& name, const Time& runTime, const dictionary& dict)
        : temporalMax(name, runTime, dict) {}
    };
    defineTypeNameAndDebug(coarseMeshMaxAlias, 0);
    addToRunTimeSelectionTable(functionObject, coarseMeshMaxAlias, dictionary);

    class coarseFieldMaxAlias : public temporalMax
    {
    public:
        TypeName("coarseFieldMax");
        coarseFieldMaxAlias(const word& name, const Time& runTime, const dictionary& dict)
        : temporalMax(name, runTime, dict) {}
    };
    defineTypeNameAndDebug(coarseFieldMaxAlias, 0);
    addToRunTimeSelectionTable(functionObject, coarseFieldMaxAlias, dictionary);

    class coarseMaxAlias : public temporalMax
    {
    public:
        TypeName("coarseMax");
        coarseMaxAlias(const word& name, const Time& runTime, const dictionary& dict)
        : temporalMax(name, runTime, dict) {}
    };
    defineTypeNameAndDebug(coarseMaxAlias, 0);
    addToRunTimeSelectionTable(functionObject, coarseMaxAlias, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::volScalarField&
Foam::functionObjects::temporalMax::lookupOrCreate
(
    const word& name,
    const dimensionSet& dims
) const
{
    if (obr_.foundObject<volScalarField>(name))
    {
        return obr_.lookupObjectRef<volScalarField>(name);
    }

    volScalarField* fieldPtr
    (
        new volScalarField
        (
            IOobject
            (
                name,
                obr_.time().name(),
                obr_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedScalar("0", dims, 0.0),
            "zeroGradient"
        )
    );
    fieldPtr->store(fieldPtr);

    return *fieldPtr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::temporalMax::temporalMax
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldNames_(),
    centre_(dict.lookupOrDefault<point>("centre", point::zero)),
    timeStart_(dict.lookupOrDefault<scalar>("timeStart", dict.lookupOrDefault<scalar>("tStart", -GREAT))),
    timeEnd_(dict.lookupOrDefault<scalar>("timeEnd", dict.lookupOrDefault<scalar>("tEnd", GREAT))),
    outputRadial_(dict.lookupOrDefault<bool>("outputRadial", true)),
    nCoarseCells_(mesh_.nCells()),
    coarseCellCentres_(mesh_.C().primitiveField()),
    cellToCoarse_(nCoarseCells_),
    coarseMax_()
{
    forAll(cellToCoarse_, i)
    {
        cellToCoarse_[i] = i;
    }

    read(dict);

    // Initialise temporal max storage for each field
    forAll(fieldNames_, i)
    {
        coarseMax_.insert(fieldNames_[i], scalarField(nCoarseCells_, -GREAT));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::temporalMax::~temporalMax()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::temporalMax::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.readIfPresent("centre", centre_);
    dict.readIfPresent("timeStart", timeStart_);
    dict.readIfPresent("tStart", timeStart_);
    dict.readIfPresent("timeEnd", timeEnd_);
    dict.readIfPresent("tEnd", timeEnd_);
    dict.readIfPresent("outputRadial", outputRadial_);

    if (dict.found("fields"))
    {
        fieldNames_ = wordList(dict.lookup("fields"));
    }
    else if (dict.found("field"))
    {
        fieldNames_ = wordList(1, dict.lookup<word>("field"));
    }
    else
    {
        fieldNames_ = wordList(1, "p");
    }

    // Ensure storage for all field names
    forAll(fieldNames_, i)
    {
        if (!coarseMax_.found(fieldNames_[i]))
        {
            coarseMax_.insert(fieldNames_[i], scalarField(nCoarseCells_, -GREAT));
        }
    }

    Info<< type() << " " << name() << ":" << nl
        << "    Tracking temporal coarse-mesh peak for fields: " << fieldNames_ << nl
        << "    Coarse cells count: " << nCoarseCells_ << nl
        << "    Time window: [" << timeStart_ << ", " << timeEnd_ << "]" << nl
        << "    Centre: " << centre_ << nl
        << endl;

    return true;
}


void Foam::functionObjects::temporalMax::topoChange
(
    const polyTopoChangeMap& map
)
{
    const labelList& cellMap = map.cellMap();
    if (cellMap.size() > 0)
    {
        labelList newCellToCoarse(mesh_.nCells(), -1);

        forAll(cellMap, i)
        {
            label oldId = cellMap[i];
            if (oldId >= 0 && oldId < cellToCoarse_.size())
            {
                newCellToCoarse[i] = cellToCoarse_[oldId];
            }
        }

        // Check reverse cell map for unrefined (merged) cells
        const labelList& rCellMap = map.reverseCellMap();
        forAll(rCellMap, i)
        {
            label index = rCellMap[i];
            if (index < -1)
            {
                label celli = -index - 2;
                if (celli >= 0 && celli < newCellToCoarse.size())
                {
                    if (i >= 0 && i < cellToCoarse_.size() && cellToCoarse_[i] >= 0)
                    {
                        newCellToCoarse[celli] = cellToCoarse_[i];
                    }
                }
            }
        }

        cellToCoarse_.transfer(newCellToCoarse);
    }
}


void Foam::functionObjects::temporalMax::mapMesh
(
    const polyMeshMap&
)
{}


bool Foam::functionObjects::temporalMax::execute()
{
    const scalar t = obr_.time().value();
    if (t < timeStart_ || t > timeEnd_)
    {
        return true;
    }

    forAll(fieldNames_, fieldi)
    {
        const word& fName = fieldNames_[fieldi];
        if (mesh_.foundObject<volScalarField>(fName))
        {
            const volScalarField& f = mesh_.lookupObject<volScalarField>(fName);
            scalarField& cMax = coarseMax_[fName];

            forAll(f, celli)
            {
                label coarseId = (celli < cellToCoarse_.size()) ? cellToCoarse_[celli] : celli;
                if (coarseId >= 0 && coarseId < nCoarseCells_)
                {
                    if (f[celli] > cMax[coarseId])
                    {
                        cMax[coarseId] = f[celli];
                    }
                }
            }
        }
    }

    return true;
}


bool Foam::functionObjects::temporalMax::write()
{
    const Time& runTime = obr_.time();
    const fileName postProcDir = runTime.globalPath()/"postProcessing"/name()/runTime.name();
    mkDir(postProcDir);

    forAll(fieldNames_, fieldi)
    {
        const word& fName = fieldNames_[fieldi];
        const scalarField& localMax = coarseMax_[fName];

        // Gather all coarse cell data across parallel processes
        List<scalarField> gatheredMax(Pstream::nProcs());
        gatheredMax[Pstream::myProcNo()] = localMax;
        Pstream::gatherList(gatheredMax);

        List<pointField> gatheredCentres(Pstream::nProcs());
        gatheredCentres[Pstream::myProcNo()] = coarseCellCentres_;
        Pstream::gatherList(gatheredCentres);

        if (Pstream::master())
        {
            scalarList globalMaxVals;
            List<point> globalCentres;
            labelList globalCoarseIDs;

            for (label procI = 0; procI < Pstream::nProcs(); ++procI)
            {
                const scalarField& pMax = gatheredMax[procI];
                const pointField& pCentres = gatheredCentres[procI];

                forAll(pMax, cI)
                {
                    globalMaxVals.append(pMax[cI]);
                    globalCentres.append(pCentres[cI]);
                    globalCoarseIDs.append(cI);
                }
            }

            // Calculate radial distance for each point
            scalarList radii(globalCentres.size());
            forAll(globalCentres, i)
            {
                radii[i] = mag(globalCentres[i] - centre_);
            }

            // Sort points by ascending radius
            labelList sortedIndices;
            sortedOrder(radii, sortedIndices);

            // Write both <field>_temporalMax.dat and <field>_coarseMax.dat for backwards compatibility
            wordList outNames = { fName + "_temporalMax.dat", fName + "_coarseMax.dat" };

            forAll(outNames, oni)
            {
                const fileName datFile = postProcDir/outNames[oni];
                OFstream os(datFile);
                if (os.good())
                {
                    os  << "# Temporal Maximum Field Output" << nl
                        << "# Field: " << fName << nl
                        << "# Time: " << runTime.name() << nl
                        << "# Time Window: [" << timeStart_ << ", " << timeEnd_ << "]" << nl
                        << "# Centre: " << centre_ << nl
                        << "# Columns:" << nl
                        << "# 1:coarseCellID  2:radius[m]  3:x[m]  4:y[m]  5:z[m]  6:"
                        << fName << "_max" << nl;

                    forAll(sortedIndices, si)
                    {
                        label idx = sortedIndices[si];
                        os  << globalCoarseIDs[idx] << " "
                            << radii[idx] << " "
                            << globalCentres[idx].x() << " "
                            << globalCentres[idx].y() << " "
                            << globalCentres[idx].z() << " "
                            << globalMaxVals[idx] << nl;
                    }

                    if (oni == 0)
                    {
                        Info<< type() << ": written fine-grained temporal peak data to " << datFile << endl;
                    }
                }
            }
        }

        // Register volScalarField on mesh
        if (mesh_.foundObject<volScalarField>(fName))
        {
            const volScalarField& f = mesh_.lookupObject<volScalarField>(fName);
            volScalarField& maxField = lookupOrCreate(fName + "TemporalMax", f.dimensions());

            forAll(maxField, celli)
            {
                label coarseId = (celli < cellToCoarse_.size()) ? cellToCoarse_[celli] : celli;
                if (coarseId >= 0 && coarseId < nCoarseCells_)
                {
                    maxField[celli] = localMax[coarseId];
                }
            }
            maxField.write();
        }
    }

    return true;
}


// ************************************************************************* //
