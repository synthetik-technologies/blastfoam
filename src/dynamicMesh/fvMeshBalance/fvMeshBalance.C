/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
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

#include "fvMeshBalance.H"
#include "decompositionMethod.H"
#include "addToRunTimeSelectionTable.H"
#include "RefineBalanceMeshObject.H"
#include "parcelCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMeshBalance, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshBalance::fvMeshBalance(fvMesh& mesh)
:
    mesh_(mesh),
    decompositionDict_
    (
        IOdictionary
        (
            IOobject
            (
                "decomposeParDict",
                mesh.time().system(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        )
    ),
    distributor_(mesh_),
    balance_(false),
    allowableImbalance_(0.2)
{
    read(decompositionDict_);
}


Foam::fvMeshBalance::fvMeshBalance
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    decompositionDict_
    (
        IOdictionary
        (
            IOobject
            (
                "decomposeParDict",
                mesh.time().system(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        )
    ),
    distributor_(mesh_),
    balance_(false),
    allowableImbalance_(0.2)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshBalance::~fvMeshBalance()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshBalance::read(const dictionary& balanceDict)
{
    if (!Pstream::parRun())
    {
        balance_ = false;
        return;
    }

    balance_ = balanceDict.lookupOrDefault("balance", true);

    if (!balance_)
    {
        return;
    }

    // Change decomposition method if entry is present
    word method
    (
        balanceDict.lookupOrDefault
        (
            "method",
            decompositionDict_.lookup<word>("method")
        )
    );
    decompositionDict_.set("method", method);

    // Add refinementHistory constraint
    if (!decompositionDict_.found("constraints"))
    {
        decompositionDict_.add("constraints", dictionary());
    }
    {
        dictionary refinementHistoryDict;
        refinementHistoryDict.add("type", "hexRefRefinementHistory");
        dictionary constraintsDict;
        decompositionDict_.subDict("constraints").set
        (
            "refinementHistory",
            refinementHistoryDict
        );

    }

    if (!decomposer_.valid())
    {
        decomposer_ = decompositionMethod::New(decompositionDict_);
    }
    if (!decomposer_->parallelAware())
    {
        FatalErrorInFunction
            << "You have selected decomposition method "
            << decomposer_->typeName
            << " which is not parallel aware." << endl
            << "Please select one that is (hierarchical, ptscotch)"
            << exit(FatalError);
    }

    if (balanceDict.found("allowableImbalance"))
    {
        allowableImbalance_ =
            balanceDict.lookup<scalar>("allowableImbalance");
    }
}


bool Foam::fvMeshBalance::canBalance() const
{
    if (!balance_)
    {
        return false;
    }

    //First determine current level of imbalance - do this for all
    // parallel runs with a changing mesh, even if balancing is disabled
    label nGlobalCells = returnReduce(mesh_.nCells(), sumOp<label>());
    scalar idealNCells =
        scalar(nGlobalCells)/scalar(Pstream::nProcs());
    scalar localImbalance = mag(scalar(mesh_.nCells()) - idealNCells);
    Foam::reduce(localImbalance, maxOp<scalar>());
    scalar maxImbalance = localImbalance/idealNCells;

    Info<<"Maximum imbalance = " << 100*maxImbalance << " %" << endl;

    //If imbalanced, construct weighted coarse graph (level 0) with node
    // weights equal to their number of subcells. This partitioning works
    // as long as the number of level 0 cells is several times greater than
    // the number of processors.
    if (maxImbalance > allowableImbalance_)
    {
        return true;
    }
    return false;
}


Foam::autoPtr<Foam::mapDistributePolyMesh>
Foam::fvMeshBalance::distribute()
{
    //Correct values on all coupled patches
    correctBoundaries<volScalarField>();
    correctBoundaries<volVectorField>();
    correctBoundaries<volSphericalTensorField>();
    correctBoundaries<volSymmTensorField>();
    correctBoundaries<volTensorField>();

    correctBoundaries<pointScalarField>();
    correctBoundaries<pointVectorField>();
    correctBoundaries<pointSphericalTensorField>();
    correctBoundaries<pointSymmTensorField>();
    correctBoundaries<pointTensorField>();

    // Decompose the mesh with uniform weights
    // The refinementHistory constraint is applied internally
    labelList finalDecomp = decomposer_().decompose
    (
        mesh_,
        scalarField(mesh_.nCells(), 1.0)
    );

    Info<< "Distributing the mesh ..." << endl;
    autoPtr<mapDistributePolyMesh> map =
        distributor_.distribute(finalDecomp);

    Info << "Successfully distributed mesh" << endl;

    scalarList procLoadNew (Pstream::nProcs(), 0.0);
    procLoadNew[Pstream::myProcNo()] = mesh_.nCells();

    reduce(procLoadNew, sumOp<List<scalar> >());

    scalar overallLoadNew = sum(procLoadNew);
    scalar averageLoadNew = overallLoadNew/double(Pstream::nProcs());

    Info << "Max deviation: " << max(Foam::mag(procLoadNew-averageLoadNew)/averageLoadNew)*100.0 << " %" << endl;

    //Correct values on all coupled patches
    correctBoundaries<volScalarField>();
    correctBoundaries<volVectorField>();
    correctBoundaries<volSphericalTensorField>();
    correctBoundaries<volSymmTensorField>();
    correctBoundaries<volTensorField>();

    correctBoundaries<pointScalarField>();
    correctBoundaries<pointVectorField>();
    correctBoundaries<pointSphericalTensorField>();
    correctBoundaries<pointSymmTensorField>();
    correctBoundaries<pointTensorField>();

    //- Update objects stored on the mesh db
    BalanceMeshObject::updateObjects(mesh_);

    //- Update objects stores on the time db
    BalanceMeshObject::updateObjects(mesh_.time());

    return map;
}

// ************************************************************************* //
