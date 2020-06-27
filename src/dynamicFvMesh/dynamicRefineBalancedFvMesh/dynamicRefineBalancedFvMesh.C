/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
     \\/     M anipulation  |
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

#include "dynamicRefineBalancedFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceInterpolate.H"
#include "volFields.H"
#include "polyTopoChange.H"
#include "surfaceFields.H"
#include "syncTools.H"
#include "pointFields.H"
#include "fvCFD.H"
#include "volPointInterpolation.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicRefineBalancedFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicRefineBalancedFvMesh, IOobject);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
Foam::label Foam::dynamicRefineBalancedFvMesh::topParentID(label p)
{
    label nextP = meshCutter()->history().splitCells()[p].parent_;
    if( nextP < 0 )
    {
        return p;
    }
    else
    {
        return topParentID(nextP);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicRefineBalancedFvMesh::dynamicRefineBalancedFvMesh
(
    const IOobject& io
)
:
    dynamicRefineFvMesh(io),
    rebalance_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicRefineBalancedFvMesh::~dynamicRefineBalancedFvMesh()
{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicRefineBalancedFvMesh::update()
{

    //Part 1 - Call normal update from dynamicRefineFvMesh
    bool hasChanged = dynamicRefineFvMesh::update();

    if( Pstream::parRun() && hasChanged)
    {
        //Correct values on all coupled patches
        correctBoundaries<scalar>();
        correctBoundaries<vector>();
        correctBoundaries<sphericalTensor>();
        correctBoundaries<symmTensor>();
        correctBoundaries<tensor>();
    }

    dictionary decomposeParDict
    (
        IOdictionary
        (
            IOobject
            (
                "decomposeParDict",
                time().system(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        )
    );
    
    rebalance_ = false; 

    // Part 2 - Load Balancing
    {    
        dictionary refineDict
        (
            IOdictionary
            (
                IOobject
                (
                    "dynamicMeshDict",
                    time().constant(),
                    *this,
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE,
                    false
                )
            ).subDict("dynamicRefineFvMeshCoeffs")
        );

        Switch enableBalancing = refineDict.lookup("enableBalancing");

        if ( Pstream::parRun() && hasChanged )
        {
            const scalar allowableImbalance =
                readScalar(refineDict.lookup("allowableImbalance"));

            //First determine current level of imbalance - do this for all
            // parallel runs with a changing mesh, even if balancing is disabled
            label nGlobalCells = globalData().nTotalCells();
            scalar idealNCells = scalar(nGlobalCells)/scalar(Pstream::nProcs());
            scalar localImbalance = mag(scalar(nCells()) - idealNCells);
            Foam::reduce(localImbalance, maxOp<scalar>());
            scalar maxImbalance = localImbalance/idealNCells;

            Info<<"Maximum imbalance = " << 100*maxImbalance << " %" << endl;

            //If imbalanced, construct weighted coarse graph (level 0) with node
            // weights equal to their number of subcells. This partitioning works
            // as long as the number of level 0 cells is several times greater than
            // the number of processors.
            if( maxImbalance > allowableImbalance && enableBalancing)
            {
                Info << "\n**Solver hold for redistribution at time = "  << time().timeName() << " s" << endl;

                rebalance_ = true;

                const labelIOList& cellLevel = meshCutter()->cellLevel();
                Map<label> coarseIDmap(100);
                labelList uniqueIndex(nCells(),0);

                label nCoarse = 0;

                forAll(cells(), cellI)
                {
                    if( cellLevel[cellI] > 0 )
                    {
                        //YO- 2D refinement uses fixed lists with unset parents;
                        //    we need to check that the parentIndex is set
                        label parentI = meshCutter()->history().parentIndex(cellI);

                        if (parentI >= 0)
                        {
                            uniqueIndex[cellI] = nCells() + topParentID
                            (
                                meshCutter()->history().parentIndex(cellI)
                            );
                        }
                        else
                        {
                            uniqueIndex[cellI] = cellI;
                        }
                        //-YO
                    }
                    else
                    {
                        uniqueIndex[cellI] = cellI;
                    }

                    if( coarseIDmap.insert(uniqueIndex[cellI], nCoarse) )
                    {
                        ++nCoarse;
                    }
                }

                // Convert to local sequential indexing and calculate coarse
                // points and weights
                labelList localIndex(nCells(),0);
                pointField coarsePoints(nCoarse,vector::zero);
                scalarField coarseWeights(nCoarse,0.0);
                label nRefinementDimensions(nGeometricD());

                forAll(uniqueIndex, cellI)
                {
                    localIndex[cellI] = coarseIDmap[uniqueIndex[cellI]];

                    // If 2D refinement (quadtree) is ever implemented, this '3'
                    // should be set in general as the number of refinement
                    // dimensions.
                    label w = (1 << (nRefinementDimensions*cellLevel[cellI]));

                    coarseWeights[localIndex[cellI]] += 1.0;
                    coarsePoints[localIndex[cellI]] += C()[cellI]/w;
                }
            
                // Set up decomposer - a separate dictionary is used here so
                // you can use a simple partitioning for decomposePar and
                // ptscotch for the rebalancing (or any chosen algorithms)
                autoPtr<decompositionMethod> decomposer
                (
                    decompositionMethod::New
                    (
                        IOdictionary
                        (
                            IOobject
                            (
                                "balanceParDict",
                                time().system(),
                                *this,
                                IOobject::MUST_READ_IF_MODIFIED,
                                IOobject::NO_WRITE
                            )
                        )
                    )
                );

                labelList finalDecomp = decomposer().decompose
                (
                    *this,
                    localIndex,
                    coarsePoints,
                    coarseWeights
                );

                scalar tolDim = globalMeshData::matchTol_ * bounds().mag();

                Info<< "Distributing the mesh ..." << endl;
                fvMeshDistribute distributor(*this, tolDim);

                Info<< "Mapping the fields ..." << endl;
                autoPtr<mapDistributePolyMesh> map =
                    distributor.distribute(finalDecomp);

                Info<< "Distribute the map ..." << endl;
                meshCutter_->distribute(map);


                Info << "Successfully distributed mesh" << endl;

                scalarList procLoadNew (Pstream::nProcs(), 0.0);
                procLoadNew[Pstream::myProcNo()] = this->nCells();

                reduce(procLoadNew, sumOp<List<scalar> >());

                scalar overallLoadNew = sum(procLoadNew);
                scalar averageLoadNew = overallLoadNew/double(Pstream::nProcs());

                Info << "New distribution: " << procLoadNew << endl;
                Info << "Max deviation: " << max(Foam::mag(procLoadNew-averageLoadNew)/averageLoadNew)*100.0 << " %" << endl;
            }
        }
    }

    if( Pstream::parRun() && rebalance_)
    {
        //Correct values on all coupled patches
        correctBoundaries<scalar>();
        correctBoundaries<vector>();
        correctBoundaries<sphericalTensor>();
        correctBoundaries<symmTensor>();
        correctBoundaries<tensor>();
    }

    return hasChanged;
}


// ************************************************************************* //
