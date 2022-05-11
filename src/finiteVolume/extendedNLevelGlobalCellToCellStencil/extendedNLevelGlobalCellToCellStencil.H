/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.


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
    Foam::extendedNLevelGlobalCellToCellStencil

Description
    Creates an n-level cell stencil using global indicies

SourceFiles
    extendedNLevelGlobalCellToCellStencil.C

\*---------------------------------------------------------------------------*/

#ifndef extendedNLevelGlobalCellToCellStencil_H
#define extendedNLevelGlobalCellToCellStencil_H

#include "MeshObject.H"
#include "PackedBoolList.H"
#include "labelPair.H"
#include "globalIndex.H"
#include "mapDistribute.H"
#include "DynamicList.H"
#include "volFields.H"
#include "cellStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class extendedNLevelGlobalCellToCellStencil Declaration
\*---------------------------------------------------------------------------*/

template<class StencilType>
class extendedNLevelGlobalCellToCellStencil
:
    public MeshObject
    <
        polyMesh,
        UpdateableMeshObject,
        extendedNLevelGlobalCellToCellStencil<StencilType>
    >
{
public:

    // Private Data

        //- Const reference to the mesh
        const polyMesh& mesh_;

        //- Number of time to extend neighbour list
        mutable label nLevels_;

        //- Global indexing
        mutable autoPtr<globalIndex> gIndexPtr_;

        //- Map of neighbours (global indexing)
        mutable Map<cellStencil> stencilMap_;

        //- List of cellCell neighbours (global)
        //  Used for distributing
        mutable labelListList cellCells_;

        mutable autoPtr<mapDistribute> mapPtr_;

        //- Calculate the cells neighbours up to the kth level
        void addCellNeighbours
        (
            const Map<labelList>& cellCells,
            const label celli,
            Map<label>& visited,
            DynamicList<label>& neighbours,
            const label level = 0
        ) const;

        void makeGlobal
        (
            const globalIndex& cfGlobalIndex,
            labelList& indices
        ) const;

        //- Combine operator for labelLists
        class unionEqOp
        {
            public:
            void operator()(labelList& x, const labelList& y) const;
        };


public:

    //- Runtime type information
    ClassName("extendedNLevelGlobalCellToCellStencil");


    // Constructors

        //- Construct from polyMesh and dictionary
        extendedNLevelGlobalCellToCellStencil
        (
            const polyMesh&,
            const label level = 1
        );

    //- Destructor
    virtual ~extendedNLevelGlobalCellToCellStencil();


    // Member Functions

        //- Build a mapDistribute for the supplied destination processor data
        static autoPtr<mapDistribute> buildMap(const List<label>& toProc);

        //- Update after refining and balancing
        virtual void updateStencil() const;

        //- Set the size of the stencil to build
        //  The largest need stencil will be build, but a smaller stencil
        // can be accessed
        virtual void setNLevel(const label n) const
        {
            nLevels_ = max(nLevels_, n);
        }

        //- Callback for geometry motion
        virtual bool movePoints();

        virtual void updateMesh(const mapPolyMesh& mpm);

        inline const globalIndex& globalNumbering() const
        {
            if (!gIndexPtr_.valid())
            {
                updateStencil();
            }
            return gIndexPtr_();
        }

        inline const Map<cellStencil>& cellCellMap() const
        {
            if (!mapPtr_.valid())
            {
                updateStencil();
            }
            return stencilMap_;
        }

        inline const cellStencil& cellCell(const label celli) const
        {
            if (!mapPtr_.valid())
            {
                updateStencil();
            }
            return stencilMap_[gIndexPtr_->toGlobal(celli)];
        }

        //- Collect data across all processors
        template<class Type>
        void collectData
        (
            const Field<Type>& fld,
            List<List<Type>>& stencilFld
        ) const;

        //- Collect data from the master cell
        template<class Type>
        void collectOwnerData
        (
            const Field<Type>& fld,
            const label level,
            Map<Type>& mapFld
        ) const;

        //- Reduce data data from the master cell
        template<class Type, class BinaryOp>
        void reduce
        (
            const Map<Type>& mapFld,
            List<Type>& fld,
            const BinaryOp& bop
        ) const;

        //- Sum surface field contributions to create cell values
        template<class Type, class WeightType>
        tmp
        <
            GeometricField
            <
                typename outerProduct<WeightType, Type>::type,
                fvPatchField,
                volMesh
            >
        > weightedSum
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            const List<List<WeightType>>& stencilWeights
        ) const;

        virtual bool writeData(Ostream&) const
        {
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "extendedNLevelGlobalCellToCellStencil.C"
    #include "extendedNLevelGlobalCellToCellStencilTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
