/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 OpenCFD Ltd.
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
    Foam::cellCellStencils::trackingInverseDistance

Description
    Inverse-distance-weighted interpolation stencil.

SourceFiles
    trackingInverseDistanceCellCellStencil.C

\*---------------------------------------------------------------------------*/

#ifndef cellCellStencils_trackingInverseDistance_H
#define cellCellStencils_trackingInverseDistance_H

#include "inverseDistanceCellCellStencil.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class voxelMeshSearch;

namespace cellCellStencils
{

/*---------------------------------------------------------------------------*\
                   Class trackingInverseDistance Declaration
\*---------------------------------------------------------------------------*/

class trackingInverseDistance
:
    public inverseDistance
{
protected:

    // Protected data

        //- Global indexing for cells
        globalIndex globalCells_;

        //- Current (global) donor cell
        labelList globalDonor_;

        //- Subset according to zone
        PtrList<fvMeshSubset> meshParts_;



    // Protected Member Functions

        //- Mark voxels of patchTypes with type of patch face
        static bool markBoundaries
        (
            const fvMesh& mesh,
            const vector& smallVec,

            const boundBox& bb,
            labelVector& nDivs,
            PackedList<2>& patchTypes,

            const labelList& cellMap,
            labelList& patchCellTypes
        );

        //- Mark all cells overlapping (a voxel covered by) a src patch
        //  with type HOLE
        void markPatchesAsHoles
        (
            PstreamBuffers& pBufs,

            // Voxel mesh for patches
            const List<treeBoundBoxList>& patchBb,
            const List<labelVector>& patchDivisions,
            const PtrList<PackedList<2>>& patchParts,

            const label srcI,
            const label tgtI,
            labelList& allCellTypes
        ) const;

        void markDonors
        (
            PstreamBuffers& pBufs,
            const List<treeBoundBoxList>& meshBb,
            const PtrList<voxelMeshSearch>& meshSearches,
            const labelList& allCellTypes,

            const label srcI,
            const label tgtI,
            labelListList& allStencil,
            labelList& allDonor
        ) const;


private:

    // Private Member Functions

        //- No copy construct
        trackingInverseDistance(const trackingInverseDistance&) = delete;

        //- No copy assignment
        void operator=(const trackingInverseDistance&) = delete;


public:

    //- Runtime type information
    TypeName("trackingInverseDistance");


    // Constructors

        //- Construct from fvMesh
        trackingInverseDistance(const fvMesh&, const dictionary&, const bool);


    //- Destructor
    virtual ~trackingInverseDistance();


    // Member Functions

        //- Update stencils. Return false if nothing changed.
        virtual bool update();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cellCellStencils
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
