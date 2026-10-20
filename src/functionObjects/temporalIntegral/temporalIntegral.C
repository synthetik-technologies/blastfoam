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

#include "temporalIntegral.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "Time.H"
#include "OFstream.H"
#include "Pstream.H"
#include "polyMeshMap.H"
#include "OSspecific.H"
#include "ListOps.H"
#include "globalIndex.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(temporalIntegral, 0);
    addToRunTimeSelectionTable(functionObject, temporalIntegral, dictionary);

    class timeIntegralAlias : public temporalIntegral
    {
    public:
        TypeName("timeIntegral");
        timeIntegralAlias(const word& name, const Time& runTime, const dictionary& dict)
        : temporalIntegral(name, runTime, dict) {}
    };
    defineTypeNameAndDebug(timeIntegralAlias, 0);
    addToRunTimeSelectionTable(functionObject, timeIntegralAlias, dictionary);

    class coarseMeshIntegralAlias : public temporalIntegral
    {
    public:
        TypeName("coarseMeshIntegral");
        coarseMeshIntegralAlias(const word& name, const Time& runTime, const dictionary& dict)
        : temporalIntegral(name, runTime, dict) {}
    };
    defineTypeNameAndDebug(coarseMeshIntegralAlias, 0);
    addToRunTimeSelectionTable(functionObject, coarseMeshIntegralAlias, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::volScalarField&
Foam::functionObjects::temporalIntegral::lookupOrCreate
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

Foam::functionObjects::temporalIntegral::temporalIntegral
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
    refValue_(dict.lookupOrDefault<scalar>("refValue", dict.lookupOrDefault<scalar>("pRef", 0.0))),
    positiveOnly_(dict.lookupOrDefault<bool>("positiveOnly", dict.found("pRef"))),
    outputRadial_(dict.lookupOrDefault<bool>("outputRadial", true)),
    nCoarseCells_(returnReduce(mesh_.nCells(), sumOp<label>())),
    coarseCellCentres_(),
    coarseCellVolumes_(),
    cellToCoarse_(mesh_.nCells()),
    coarseIntegral_()
{
    const globalIndex coarseIndex(mesh_.nCells());
    forAll(cellToCoarse_, i)
    {
        cellToCoarse_[i] = coarseIndex.toGlobal(i);
    }

    List<pointField> centres(Pstream::nProcs());
    centres[Pstream::myProcNo()] = mesh_.C().primitiveField();
    Pstream::gatherList(centres);
    if (Pstream::master())
    {
        coarseCellCentres_ =
            ListListOps::combine<pointField>(centres, accessOp<pointField>());
    }

    List<scalarField> volumes(Pstream::nProcs());
    volumes[Pstream::myProcNo()] = mesh_.V().primitiveField();
    Pstream::gatherList(volumes);
    Pstream::scatterList(volumes);
    coarseCellVolumes_ =
        ListListOps::combine<scalarField>(volumes, accessOp<scalarField>());

    read(dict);

    // Initialise temporal integral storage for each field
    forAll(fieldNames_, i)
    {
        coarseIntegral_.insert(fieldNames_[i], scalarField(nCoarseCells_, 0.0));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::temporalIntegral::~temporalIntegral()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::temporalIntegral::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.readIfPresent("centre", centre_);
    dict.readIfPresent("timeStart", timeStart_);
    dict.readIfPresent("tStart", timeStart_);
    dict.readIfPresent("timeEnd", timeEnd_);
    dict.readIfPresent("tEnd", timeEnd_);
    dict.readIfPresent("outputRadial", outputRadial_);
    dict.readIfPresent("refValue", refValue_);
    dict.readIfPresent("pRef", refValue_);
    dict.readIfPresent("positiveOnly", positiveOnly_);

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
        if (!coarseIntegral_.found(fieldNames_[i]))
        {
            coarseIntegral_.insert(fieldNames_[i], scalarField(nCoarseCells_, 0.0));
        }
    }

    Info<< type() << " " << name() << ":" << nl
        << "    Tracking temporal integral for fields: " << fieldNames_ << nl
        << "    Reference value: " << refValue_ << nl
        << "    Positive only: " << positiveOnly_ << nl
        << "    Coarse cells count: " << nCoarseCells_ << nl
        << "    Time window: [" << timeStart_ << ", " << timeEnd_ << "]" << nl
        << "    Centre: " << centre_ << nl
        << endl;

    return true;
}


void Foam::functionObjects::temporalIntegral::topoChange
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


void Foam::functionObjects::temporalIntegral::mapMesh
(
    const polyMeshMap&
)
{}


void Foam::functionObjects::temporalIntegral::distribute
(
    const polyDistributionMap& map
)
{
    map.distributeCellData(cellToCoarse_);
}


bool Foam::functionObjects::temporalIntegral::execute()
{
    const scalar t = obr_.time().value();
    if (t < timeStart_ || t > timeEnd_)
    {
        return true;
    }

    const scalar dt = obr_.time().deltaTValue();
    if (dt <= 0)
    {
        return true;
    }

    const scalarField& V = mesh_.V().primitiveField();

    forAll(fieldNames_, fieldi)
    {
        const word& fName = fieldNames_[fieldi];
        if (mesh_.foundObject<volScalarField>(fName))
        {
            const volScalarField& f = mesh_.lookupObject<volScalarField>(fName);
            scalarField& cInt = coarseIntegral_[fName];

            // Volume average over the coarse cell (refinement keeps its
            // volume), summed over the processors on write
            forAll(f, celli)
            {
                label coarseId = (celli < cellToCoarse_.size()) ? cellToCoarse_[celli] : celli;
                if (coarseId >= 0 && coarseId < nCoarseCells_)
                {
                    scalar val = f[celli] - refValue_;
                    if (positiveOnly_)
                    {
                        val = max(val, 0.0);
                    }
                    cInt[coarseId] +=
                        val*V[celli]/coarseCellVolumes_[coarseId]*dt;
                }
            }
        }
    }

    return true;
}


bool Foam::functionObjects::temporalIntegral::write()
{
    const Time& runTime = obr_.time();
    const fileName postProcDir = runTime.globalPath()/"postProcessing"/name()/runTime.name();
    mkDir(postProcDir);

    forAll(fieldNames_, fieldi)
    {
        const word& fName = fieldNames_[fieldi];
        // Sum over the processors (the fine cells of a coarse cell can be
        // on several processors after redistribution)
        scalarField localInt(coarseIntegral_[fName]);
        Pstream::listCombineGather(localInt, plusEqOp<scalar>());
        Pstream::listCombineScatter(localInt);

        if (Pstream::master())
        {
            const scalarField& globalIntVals = localInt;
            const pointField& globalCentres = coarseCellCentres_;
            const labelList globalCoarseIDs(identityMap(nCoarseCells_));

            // Calculate radial distance for each point
            scalarList radii(globalCentres.size());
            forAll(globalCentres, i)
            {
                radii[i] = mag(globalCentres[i] - centre_);
            }

            // Sort points by ascending radius
            labelList sortedIndices;
            sortedOrder(radii, sortedIndices);

            const fileName datFile = postProcDir/(fName + "_temporalIntegral.dat");
            OFstream os(datFile);
            if (os.good())
            {
                os  << "# Temporal Integral Field Output" << nl
                    << "# Field: " << fName << nl
                    << "# Time: " << runTime.name() << nl
                    << "# Time Window: [" << timeStart_ << ", " << timeEnd_ << "]" << nl
                    << "# Centre: " << centre_ << nl
                    << "# Reference Value: " << refValue_ << nl
                    << "# Positive Only: " << positiveOnly_ << nl
                    << "# Columns:" << nl
                    << "# 1:coarseCellID  2:radius[m]  3:x[m]  4:y[m]  5:z[m]  6:"
                    << fName << "_integral" << nl;

                forAll(sortedIndices, si)
                {
                    label idx = sortedIndices[si];
                    os  << globalCoarseIDs[idx] << " "
                        << radii[idx] << " "
                        << globalCentres[idx].x() << " "
                        << globalCentres[idx].y() << " "
                        << globalCentres[idx].z() << " "
                        << globalIntVals[idx] << nl;
                }

                Info<< type() << ": written fine-grained temporal integral data to " << datFile << endl;
            }
        }

        // Register volScalarField on mesh
        if (mesh_.foundObject<volScalarField>(fName))
        {
            const volScalarField& f = mesh_.lookupObject<volScalarField>(fName);
            volScalarField& intField = lookupOrCreate(fName + "TemporalIntegral", f.dimensions()*dimTime);

            forAll(intField, celli)
            {
                label coarseId = (celli < cellToCoarse_.size()) ? cellToCoarse_[celli] : celli;
                if (coarseId >= 0 && coarseId < nCoarseCells_)
                {
                    intField[celli] = localInt[coarseId];
                }
            }
            intField.write();
        }
    }

    return true;
}


// ************************************************************************* //
