/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

Class
    Foam::cellCellStencils::inverseDistance

Description
    Inverse-distance-weighted interpolation stencil.

    hole finding:
    - mark boundary faces on helper (voxel) mesh
    - mark any cell overlaying these voxels
    - use flood filling to find any unreachable cell
    Alternative is to use an octree of the boundary faces and determine
    directly for all cells whether we are outside. Might be slow though.

SourceFiles
    inverseDistanceCellCellStencil.C

\*---------------------------------------------------------------------------*/

#ifndef cellCellStencils_inverseDistance_H
#define cellCellStencils_inverseDistance_H

#include "cellCellStencil.H"
#include "volFields.H"
#include "labelVector.H"
#include "treeBoundBoxList.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMeshSubset;

namespace cellCellStencils
{

/*---------------------------------------------------------------------------*\
                     Class inverseDistance Declaration
\*---------------------------------------------------------------------------*/

class inverseDistance
:
    public cellCellStencil
{
protected:

    // Protected data

        //- Dictionary of motion control parameters
        const dictionary dict_;

        //- Small perturbation vector for geometric tests
        vector smallVec_;

        //- Per cell the cell type
        labelList cellTypes_;

        //- Indices of interpolated cells
        labelList interpolationCells_;

        //- Fetch interpolated cells
        autoPtr<mapDistribute> cellInterpolationMap_;

        //- Interpolation stencil
        labelListList cellStencil_;

        //- Interpolation weights
        scalarListList cellInterpolationWeights_;

        //- Amount of interpolation
        volScalarField cellInterpolationWeight_;


   // Protected Member Functions

        // Voxel representation

            //- Convert ijk indices into single index
            static label index(const labelVector& nDivs, const labelVector&);

            //- Convert single index into ijk
            static labelVector index3(const labelVector& nDivs, const label);

            //- Convert coordinate into ijk
            static labelVector index3
            (
                const boundBox& bb,
                const labelVector& nDivs,
                const point& pt
            );

            //- Convert index back into coordinate
            static point position
            (
                const boundBox& bb,
                const labelVector& nDivs,
                const label boxI
            );

            //- Fill all elements overlapping subBb with value val
            static void fill
            (
                PackedList<2>& elems,
                const boundBox& bb,
                const labelVector& nDivs,
                const boundBox& subBb,
                const unsigned int val
            );

            //- Is any voxel inside subBb set to val
            static bool overlaps
            (
                const boundBox& bb,
                const labelVector& nDivs,
                const PackedList<2>& voxels,
                const treeBoundBox& subBb,
                const unsigned int val
            );

        //- Mark voxels of patchTypes with type of patch face
        static void markBoundaries
        (
            const fvMesh& mesh,
            const vector& smallVec,

            const boundBox& bb,
            const labelVector& nDivs,
            PackedList<2>& patchTypes,
            const labelList& cellMap,
            labelList& patchCellTypes
        );

        //- Calculate bounding box of cell
        static treeBoundBox cellBb
        (
            const primitiveMesh& mesh,
            const label celli
        );

        //- Mark all cells overlapping (a voxel covered by) a src patch
        //  with type HOLE
        void markPatchesAsHoles
        (
            PstreamBuffers& pBufs,

            // Mesh bb and data
            const PtrList<fvMeshSubset>& meshParts,

            // Helper mesh for patches
            const List<treeBoundBoxList>& patchBb,
            const List<labelVector>& patchDivisions,
            const PtrList<PackedList<2>>& patchParts,

            const label srcI,
            const label tgtI,
            labelList& allCellTypes
        ) const;

        //- If multiple donors meshes: decide which is best
        bool betterDonor
        (
            const label destMesh,
            const label currentDonorMesh,
            const label newDonorMesh
        ) const;

        //- Determine donors for all tgt cells
        void markDonors
        (
            const globalIndex& globalCells,
            PstreamBuffers& pBufs,
            const PtrList<fvMeshSubset>& meshParts,
            const List<treeBoundBoxList>& meshBb,
            const labelList& allCellTypes,

            const label srcI,
            const label tgtI,
            labelListList& allStencil,
            labelList& allDonor
        ) const;

        //- Replacement of regionSplit
        void uncompactedRegionSplit
        (
            const fvMesh& mesh,
            const globalIndex& globalFaces,
            const label nZones,
            const labelList& zoneID,
            const labelList& cellTypes,
            const boolList& isBlockedFace,
            labelList& cellRegion
        ) const;
        autoPtr<globalIndex> compactedRegionSplit
        (
            const fvMesh& mesh,
            const globalIndex& globalFaces,
            labelList& cellRegion
        ) const;

        //- Do flood filling to detect unreachable (from patches) sections
        //  of mesh
        void findHoles
        (
            const globalIndex& globalCells,
            const fvMesh& mesh,
            const labelList& zoneID,
            const labelListList& stencil,
            labelList& cellTypes
        ) const;

        //- Seed faces of cell with wantedFraction (if higher than current)
        void seedCell
        (
            const label cellI,
            const scalar wantedFraction,
            PackedBoolList& isFront,
            scalarField& fraction
        ) const;

        //- Surround holes with layer(s) of interpolated cells
        void walkFront
        (
            const scalar layerRelax,
            const labelListList& allStencil,
            labelList& allCellTypes,
            scalarField& allWeight
        ) const;

        //- Create stencil starting from the donor containing the acceptor
        virtual void createStencil(const globalIndex&);


private:

    // Private Member Functions

        //- No copy construct
        inverseDistance(const inverseDistance&) = delete;

        //- No copy assignment
        void operator=(const inverseDistance&) = delete;


public:

    //- Runtime type information
    TypeName("inverseDistance");


    // Constructors

        //- Construct from fvMesh
        inverseDistance(const fvMesh&, const dictionary&, const bool);


    //- Destructor
    virtual ~inverseDistance();


    // Member Functions

        //- Update stencils. Return false if nothing changed.
        virtual bool update();

        //- Callback for geometry motion
        virtual void updateMesh(const mapPolyMesh& map);

        //- Return the cell type list
        virtual const labelUList& cellTypes() const
        {
            return cellTypes_;
        }

        //- Indices of interpolated cells
        virtual const labelUList& interpolationCells() const
        {
            return interpolationCells_;
        }

        //- Return a communication schedule
        virtual const mapDistribute& cellInterpolationMap() const
        {
            if (!cellInterpolationMap_.valid())
            {
                const_cast<inverseDistance&>(*this).update();
            }
            return cellInterpolationMap_();
        }

        //- Per interpolated cell the neighbour cells (in terms of slots as
        //  constructed by above cellInterpolationMap) to interpolate
        virtual const labelListList& cellStencil() const
        {
            return cellStencil_;
        }

        //- Weights for cellStencil
        virtual const scalarListList& cellInterpolationWeights() const
        {
            return cellInterpolationWeights_;
        }

        //- Per interpolated cell the interpolation factor. (0 = use
        //  calculated, 1 = use interpolated)
        virtual const scalarList& cellInterpolationWeight() const
        {
            return cellInterpolationWeight_;
        }

        //- Calculate inverse distance weights for a single acceptor
        virtual void stencilWeights
        (
            const point& sample,
            const pointList& donorCcs,
            scalarList& weights
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cellCellStencils
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
