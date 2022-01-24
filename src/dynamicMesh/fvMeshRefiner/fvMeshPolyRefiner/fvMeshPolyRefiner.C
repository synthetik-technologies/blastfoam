/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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

#include "fvMeshPolyRefiner.H"
#include "RefineBalanceMeshObject.H"
#include "polyTopoChange.H"
#include "parcelCloud.H"
#include "prismatic2DRefinement.H"
#include "polyhedralRefinement.H"
#include "polyRefinementConstraint.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMeshPolyRefiner, 0);
    addToRunTimeSelectionTable(fvMeshRefiner, fvMeshPolyRefiner, fvMesh);
    addToRunTimeSelectionTable(fvMeshRefiner, fvMeshPolyRefiner, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::fvMeshPolyRefiner::selectUnrefinePoints
(
    const scalar unrefineLevel,
    const PackedBoolList& markedCell,
    const scalarField& pFld
) const
{
    DynamicList<label> candidates(mesh_.nPoints());
    forAll(mesh_.points(), pointI)
    {
        if (pFld[pointI] < unrefineLevel)
        {
            // Check that all cells are not marked
            const labelList& pCells = mesh_.pointCells()[pointI];

            bool hasMarked = false;

            forAll(pCells, pCelli)
            {
                if (markedCell.get(pCells[pCelli]))
                {
                    hasMarked = true;
                    break;
                }
            }

            if (!hasMarked)
            {
                candidates.append(pointI);
            }
        }
    }

    // Guarantee 2:1 refinement after refinement
    labelList consistentSet
    (
        refiner_->consistentUnrefinement
        (
            candidates.shrink(),
            false               // Add to set to guarantee 2:1
        )
    );

    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " points to split " << mesh_.globalData().nTotalPoints()
        << "." << endl;

    return consistentSet;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshPolyRefiner::fvMeshPolyRefiner(fvMesh& mesh)
:
    fvMeshRefiner(mesh),

    refiner_(nullptr),
    isRefining_(false),
    isUnrefining_(false)
{
    // Added refinement history decomposition constraint to keep all
    // cells with the same parent together
    {
        dictionary refinementHistoryDict("refinementHistory");
        refinementHistoryDict.add
        (
            "type",
            polyRefinementConstraint::typeName
        );
        balancer_.addConstraint(refinementHistoryDict);
    }

    // Get number of valid geometric dimensions
    const label nGeometricDirs = mesh_.nGeometricD();

    switch(nGeometricDirs)
    {
        case 3:
            // Add the polyhedralRefinement engine for
            // 3D isotropic refinement
            Info<< "3D case detected. "
                << "Adding polyhedralRefinement topology modifier" << endl;
            refiner_.set
            (
                new polyhedralRefinement
                (
                    mesh,
                    dict_
                )
            );
            break;

        case 2:
            // Add the prismatic2DRefinement engine for
            // 2D isotropic refinement
            Info<< "2D case detected. "
                << "Adding prismatic2DRefinement topology modifier" << endl;
            refiner_.set
            (
                new prismatic2DRefinement
                (
                    mesh,
                    dict_
                )
            );
            break;

        case 1:
            FatalErrorInFunction
                << "1D case detected. No valid refinement strategy is"
                <<  " available for 1D cases."
                << abort(FatalError);
            break;

        default:
            FatalErrorInFunction
                << "Invalid number of geometric meshes detected: "
                << nGeometricDirs
                << nl << "It appears that this mesh is neither 1D, 2D or 3D."
                << abort(FatalError);

    }
}


Foam::fvMeshPolyRefiner::fvMeshPolyRefiner
(
    fvMesh& mesh,
    const dictionary& dict,
    const bool force,
    const bool read
)
:
    fvMeshRefiner(mesh, dict, force, read),

    refiner_(nullptr),
    isRefining_(false),
    isUnrefining_(false)
{
    // Added refinement history decomposition constraint to keep all
    // cells with the same parent together
    {
        dictionary refinementHistoryDict("refinementHistory");
        refinementHistoryDict.add
        (
            "type",
            polyRefinementConstraint::typeName
        );
        balancer_.addConstraint(refinementHistoryDict);
    }

    readDict(dict);

    // Get number of valid geometric dimensions
    const label nGeometricDirs = mesh_.nGeometricD();

    switch(nGeometricDirs)
    {
        case 3:
            // Add the polyhedralRefinement engine for
            // 3D isotropic refinement
            Info<< "3D case detected. "
                << "Adding polyhedralRefinement topology modifier" << endl;
            refiner_.set
            (
                new polyhedralRefinement
                (
                    mesh,
                    dict_,
                    read
                )
            );
            break;

        case 2:
            // Add the prismatic2DRefinement engine for
            // 2D isotropic refinement
            Info<< "2D case detected. "
                << "Adding prismatic2DRefinement topology modifier" << endl;
            refiner_.set
            (
                new prismatic2DRefinement
                (
                    mesh,
                    dict_,
                    read
                )
            );
            break;

        case 1:
            FatalErrorInFunction
                << "1D case detected. No valid refinement strategy is"
                <<  " available for 1D cases."
                << abort(FatalError);
            break;

        default:
            FatalErrorInFunction
                << "Invalid number of geometric meshes detected: "
                << nGeometricDirs
                << nl << "It appears that this mesh is neither 1D, 2D or 3D."
                << abort(FatalError);

    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshPolyRefiner::~fvMeshPolyRefiner()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMeshPolyRefiner::readDict(const dictionary& dict)
{
    fvMeshRefiner::readDict(dict);

    // Check whether the number of Refinement buffer layers is smaller than
    // 3 for 2D cases
    static bool hasWarnedRefine = false;
    if
    (
        nRefinementBufferLayers_ < 3
     && !hasWarnedRefine
     && mesh_.nGeometricD() == 2
    )
    {
        hasWarnedRefine = true;
        WarningInFunction
            << "Using " << nRefinementBufferLayers_
            << " refinement buffer layers" << nl
            << "Make sure that the number of refinement buffer layers is "
            << "at least 3 in order to avoid problems with edge level "
            << "in 2 dimensional cases"
            << endl;
        nRefinementBufferLayers_ = 3;
    }

    // Check whether the number of unrefinement buffer layers is smaller than
    // number of refinement buffer layers + 2
    static bool hasWarnedUnrefine = false;
    if (nUnrefinementBufferLayers_ < 2 && !hasWarnedUnrefine)
    {
        hasWarnedUnrefine = true;
        WarningInFunction
            << "Using " << nUnrefinementBufferLayers_
            << " unrefinement buffer layers" << nl
            << "Make sure that the number of unrefinement buffer layers is "
            << "at least 2 in order to avoid problems with edge level "
            << "inconsistency when refinement and unrefinement are performed in "
            << "same iteration."
            << endl;
        nUnrefinementBufferLayers_ = 2;
    }
}


bool Foam::fvMeshPolyRefiner::refine
(
    const scalarField& error,
    const labelList& maxCellLevel,
    const scalar lowerRefineLevel,
    const scalar upperRefineLevel,
    const scalar unrefineLevel
)
{
    bool hasChanged = false;

    if (preUpdate())
    {
        // Cells marked for refinement or otherwise protected from unrefinement.
        PackedBoolList refineCell(mesh_.nCells());

        if (canRefine(true))
        {
            labelList maxRefinement(maxCellLevel);
            setMaxCellLevel(maxRefinement);

            // Determine candidates for refinement (looking at field only)
            selectRefineCandidates
            (
                lowerRefineLevel,
                upperRefineLevel,
                error,
                refineCell
            );

            // Extend with a buffer layer to refine neighbour cells
            for (label i = 0; i < nRefinementBufferLayers_; i++)
            {
                extendMarkedCells(refineCell, maxRefinement, i == 0, force_);
            }

            forAll(protectedPatches_, patchi)
            {
                const polyPatch& p =
                    mesh_.boundaryMesh()[protectedPatches_[patchi]];
                forAll(p.faceCells(), facei)
                {
                    label own = mesh_.faceOwner()[facei + p.start()];
                    refineCell.set(own, false);
                }
            }

            // Select subset of candidates. Take into account max allowable
            // cells, refinement level, protected cells.
            labelList cellsToRefine
            (
                refiner_->consistentRefinement
                (
                    selectRefineCells
                    (
                        maxCells_,
                        maxRefinement,
                        refineCell
                    ),
                    true
                )
            );

            label nCellsToRefine = returnReduce
            (
                cellsToRefine.size(), sumOp<label>()
            );

            Info<< "Selected "
                << nCellsToRefine
                << " cells for refinement out of "
                << mesh_.globalData().nTotalCells()
                << "." << endl;

            if (nCellsToRefine > 0)
            {
                isRefining_ = true;
                autoPtr<mapPolyMesh> map = refiner_->refine
                (
                    mesh_,
                    cellsToRefine
                );

                // Update refineCell. Note that some of the marked ones have
                // not been refined due to constraints.
                {
                    const labelList& cellMap = map().cellMap();
                    const labelList& reverseCellMap = map().reverseCellMap();

                    // Create new refineCell
                    PackedBoolList newRefineCell(cellMap.size());

                    forAll(cellMap, celli)
                    {
                        label oldCelli = cellMap[celli];

                        if (oldCelli < 0)
                        {
                            newRefineCell.set(celli, 1);
                        }
                        else if (reverseCellMap[oldCelli] != celli)
                        {
                            newRefineCell.set(celli, 1);
                        }
                        else
                        {
                            newRefineCell.set(celli, refineCell.get(oldCelli));
                        }
                    }
                    refineCell.transfer(newRefineCell);
                }

                hasChanged = true;
                isRefining_ = false;
            }
        }

        if (canUnrefine(true))
        {
            // Extend with a buffer layer to not unrefine neighbour cells
            for (label i = 0; i < nUnrefinementBufferLayers_; i++)
            {
                extendMarkedCellsAcrossPoints(refineCell);
            }

            forAll(protectedPatches_, patchi)
            {
                const polyPatch& p =
                    mesh_.boundaryMesh()[protectedPatches_[patchi]];
                forAll(p.faceCells(), facei)
                {
                    label own = mesh_.faceOwner()[facei + p.start()];
                    refineCell.set(own, true);
                }
            }

            // Select unrefineable points that are not marked in refineCell
            labelList pointsToUnrefine
            (
                selectUnrefinePoints
                (
                    unrefineLevel,
                    refineCell,
                    maxCellField(error)
                )
            );

            label nSplitPoints = returnReduce
            (
                pointsToUnrefine.size(),
                sumOp<label>()
            );

            if (nSplitPoints > 0)
            {
                isUnrefining_ = true;
                hasChanged = hasChanged || refiner_->unrefine
                (
                    mesh_,
                    pointsToUnrefine
                );

                isUnrefining_ = false;
            }
        }

        reduce(hasChanged, orOp<bool>());
        if (balance())
        {
            hasChanged = true;
        }
        else if (hasChanged)
        {
            //- Update objects stored on the mesh db
            RefineMeshObject::updateObjects(mesh_);
        }
        mesh_.topoChanging(hasChanged);

        if (hasChanged)
        {
            // Reset moving flag (if any). If not using inflation we'll not
            // move, if are using inflation any follow on movePoints will set
            // it.
            mesh_.moving(false);
            mesh_.setInstance(mesh_.time().timeName());
        }
    }

    return hasChanged;
}

void Foam::fvMeshPolyRefiner::updateMesh(const mapPolyMesh& map)
{
    fvMeshRefiner::updateMesh(map);
    if (!isBalancing_)
    {
        refiner_->updateMesh(map);
    }
}


void Foam::fvMeshPolyRefiner::distribute(const mapDistributePolyMesh& map)
{
    fvMeshRefiner::distribute(map);
    refiner_->distribute(map);
}


bool Foam::fvMeshPolyRefiner::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    if (debug)
    {
        volScalarField clusters
        (
            volScalarField::New
            (
                "clusters",
                mesh_,
                dimensionedScalar(dimless, 0)
            )
        );
        labelList cellClusters;
        refiner_->getCellClusters(cellClusters);
        forAll(clusters, celli)
        {
            clusters[celli] = cellClusters[celli];
        }
        clusters.write();
    }
    return
        fvMeshRefiner::writeObject(fmt, ver, cmp, write)
     && refiner_->write();
}


// ************************************************************************* //
