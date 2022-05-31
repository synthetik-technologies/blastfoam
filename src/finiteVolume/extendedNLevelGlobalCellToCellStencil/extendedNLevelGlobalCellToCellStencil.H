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

template<class Type>
using PairMap = HashTable<Type, labelPair, Hash<labelPair>>;

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
        label nLevels_;

        //- Minimum number of cell neighbours
        labelList nNbrs_;

        //- Global indexing
        mutable autoPtr<globalIndex> gIndexPtr_;

        //- Map of neighbours (global indexing)
        mutable Map<cellStencil> stencilMap_;

        //- List of cellCell neighbours (global)
        //  Used for distributing
        mutable labelListList cellCells_;

        mutable autoPtr<mapDistribute> mapPtr_;

        //- Non local stencil cells
        mutable Map<label> nonlocalCells_;

        //- Non local stencil owners
        mutable Map<label> nonlocalOwners_;

        //- Calculate the cells neighbours up to the kth level
        void addCellNeighbours
        (
            const Map<labelList>& cellCells,
            const label celli,
            const label level,
            const label maxLevel,
            labelHashSet& neighbours,
            labelHashSet& added,
            DynamicList<label>& nbrs
        ) const;

        //- Calculate the cells neighbours up to the kth level
        labelList addCellNeighbours
        (
            const Map<labelList>& cellCells,
            const label celli,
            labelHashSet& neighbours,
            DynamicList<label>& nbrs
        ) const;

        void makeGlobal
        (
            const globalIndex& cfGlobalIndex,
            labelList& indices
        ) const;

        void createMap() const;
        const mapDistribute& map() const;

public:

    //- Runtime type information
    ClassName("extendedNLevelGlobalCellToCellStencil");


    // Constructors

        //- Construct from polyMesh and dictionary
        extendedNLevelGlobalCellToCellStencil
        (
            const polyMesh&,
            const label level
        );

        //- Construct from polyMesh and dictionary
        extendedNLevelGlobalCellToCellStencil
        (
            const polyMesh&,
            const labelList& nNbrs
        );

    //- Destructor
    virtual ~extendedNLevelGlobalCellToCellStencil();


    // Member Functions

        //- Build a mapDistribute for the supplied destination processor data
        static autoPtr<mapDistribute> buildMap(const List<label>& toProc);

        //- Update after refining and balancing
        virtual void calcStencil() const;

        //- Update local/nonlocal indices and send
        void update() const;

        //- Set the size of the stencil to build
        //  The largest need stencil will be build, but a smaller stencil
        // can be accessed
        virtual void setNLevel(const label n)
        {
            nLevels_ = max(nLevels_, n);
            nNbrs_.clear();
        }

        virtual void setNNbrs(const labelList& nNbrs)
        {
            nNbrs_ = nNbrs;
            nLevels_ = -1;
        }


        //- Callback for geometry motion
        virtual bool movePoints();

        virtual void updateMesh(const mapPolyMesh& mpm);

        inline const globalIndex& gIndex() const
        {
            if (!gIndexPtr_.valid())
            {
                calcStencil();
            }
            return gIndexPtr_();
        }

        inline const Map<cellStencil>& cellCellMap() const
        {
            if (!gIndexPtr_.valid())
            {
                calcStencil();
            }
            return stencilMap_;
        }

        inline Map<cellStencil>& cellCellMapRef()
        {
            mapPtr_.clear();
            if (!gIndexPtr_.valid())
            {
                calcStencil();
            }
            return stencilMap_;
        }

        inline const cellStencil& cellCell(const label celli) const
        {
            if (!gIndexPtr_.valid())
            {
                calcStencil();
            }
            return stencilMap_[gIndexPtr_->toGlobal(celli)];
        }

        //- Collect data across all processors
        template<class Type, template<class> class ListType>
        void collectData
        (
            const UList<Type>& fld,
            List<ListType<Type> >& stencilFld
        ) const;

        //- Collect data across all processors
        template<class Type, template<class> class ListType>
        void collectData
        (
            const PairMap<Type>& mapFld,
            List<ListType<Type> >& stencilFld
        ) const;

        //- Collect data from the master cell
        template<class Type>
        void collectOwnerData(const UList<Type>& fld, Map<Type>& mapFld) const;

        //- Collect data from the master cell
        template<class Type>
        void collectOwnerData(Map<Type>& mapFld) const;

        //- Collect data from the master cell
        template<class Type>
        void collectNbrData(Map<Type>& mapFld) const;

        //- Reduce data data to the master cell
        template<class Type, class BinaryOp>
        void reduce
        (
            const Map<Type>& mapFld,
            UList<Type>& fld,
            const BinaryOp& bop
        ) const;

        //- Sum surface field contributions to create cell values
        template<class Type, class WeightType, template<class> class ListType>
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
            const UList<ListType<WeightType>>& stencilWeights
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
