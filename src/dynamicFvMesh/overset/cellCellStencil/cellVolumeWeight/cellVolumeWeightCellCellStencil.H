/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2018 OpenCFD Ltd.
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
    Foam::cellCellStencils::cellVolumeWeight

Description
    Volume-weighted interpolation stencil

SourceFiles
    cellVolumeWeightCellCellStencil.C
    cellVolumeWeightCellCellStencilTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef cellCellStencils_cellVolumeWeight_H
#define cellCellStencils_cellVolumeWeight_H

#include "cellCellStencil.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cellCellStencils
{

/*---------------------------------------------------------------------------*\
                     Class cellVolumeWeight Declaration
\*---------------------------------------------------------------------------*/

class cellVolumeWeight
:
    public cellCellStencil
{
protected:

    // Static data members

        //- Default overlap tolerance. Fraction of volume
        static scalar defaultOverlapTolerance_;


    // Protected data

        //- Dictionary of motion control parameters
        const dictionary dict_;

        //- Tolerance for volume overlap. Fraction of volume
        scalar overlapTolerance_;

        //- Per cell the cell type
        mutable labelList cellTypes_;

        //- Indices of interpolated cells
        mutable labelList interpolationCells_;

        //- Fetch interpolated cells
        mutable autoPtr<mapDistribute> cellInterpolationMap_;

        //- Interpolation stencil
        mutable labelListList cellStencil_;

        //- Interpolation weights
        mutable List<scalarList> cellInterpolationWeights_;

        //- Amount of interpolation
        mutable volScalarField cellInterpolationWeight_;


   // Protected Member Functions

        void walkFront
        (
            const scalar layerRelax,
            labelList& allCellTypes,
            scalarField& allWeight
        ) const;

        //- Find cells next to cells of type PATCH
        void findHoles
        (
            const globalIndex& globalCells,
            const fvMesh& mesh,
            const labelList& zoneID,
            const labelListList& stencil,
            labelList& cellTypes
        ) const;

        //- according to additionalDocumentation/MEJ_oversetMesh.txt
        void markPatchCells
        (
            const fvMesh& mesh,
            const labelList& cellMap,
            labelList& patchCellTypes
        ) const;

        void combineCellTypes
        (
            const label subZoneID,
            const fvMesh& subMesh,
            const labelList& subCellMap,

            const label donorZoneID,
            const labelListList& toOtherCells,
            const List<scalarList>& weights,
            const labelList& otherCells,
            const labelList& interpolatedOtherPatchTypes,

            labelListList& allStencil,
            scalarListList& allWeights,
            labelList& allCellTypes,
            labelList& allDonorID
        ) const;

        //- interpolate (= combine) patch types
        void interpolatePatchTypes
        (
            const labelListList& addressing,
            const labelList& patchTypes,
            labelList& result
        ) const;

        //- interpolate (= combine) patch types
        void interpolatePatchTypes
        (
            const autoPtr<mapDistribute>& mapPtr,
            const labelListList& addressing,
            const labelList& patchTypes,
            labelList& result
        ) const;


private:

    // Private Member Functions

        //- No copy construct
        cellVolumeWeight(const cellVolumeWeight&) = delete;

        //- No copy assignment
        void operator=(const cellVolumeWeight&) = delete;


public:

    //- Runtime type information
    TypeName("cellVolumeWeight");


    // Constructors

        //- Construct from fvMesh
        cellVolumeWeight(const fvMesh&, const dictionary&);


    //- Destructor
    virtual ~cellVolumeWeight();


    // Member Functions

        //- Access to volume overlap tolerance
        scalar overlapTolerance() const
        {
            return overlapTolerance_;
        }

        //- Update stencils. Return false if nothing changed.
        virtual bool update() const;

        //- Callback for geometry motion
        virtual void updateMesh(const mapPolyMesh& map);

        //- Return the cell type list
        virtual const labelUList& cellTypes() const
        {
            if (needUpdate_)
            {
                update();
            }
            return cellTypes_;
        }

        //- Indices of interpolated cells
        virtual const labelUList& interpolationCells() const
        {
            if (needUpdate_)
            {
                update();
            }
            return interpolationCells_;
        }

        //- Return a communication schedule
        virtual const mapDistribute& cellInterpolationMap() const
        {
            if (!cellInterpolationMap_.valid() || needUpdate_)
            {
                update();
            }
            return cellInterpolationMap_();
        }

        //- Per interpolated cell the neighbour cells (in terms of slots as
        //  constructed by above cellInterpolationMap) to interpolate
        virtual const labelListList& cellStencil() const
        {
            if (needUpdate_)
            {
                update();
            }
            return cellStencil_;
        }

        //- Weights for cellStencil
        virtual const List<scalarList>& cellInterpolationWeights() const
        {
            if (needUpdate_)
            {
                update();
            }
            return cellInterpolationWeights_;
        }

        //- Per interpolated cell the interpolation factor. (0 = use
        //  calculated, 1 = use interpolated)
        virtual const scalarList& cellInterpolationWeight() const
        {
            if (needUpdate_)
            {
                update();
            }
            return cellInterpolationWeight_;
        }

        //- Calculate inverse distance weights for a single acceptor. Revert
        //  to inverse distance (so not consistent with volume overlap!)
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
