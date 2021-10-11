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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class cellStencil;
Ostream& operator<<(Ostream&, const cellStencil&);
Istream& operator>>(Istream&, cellStencil&);

class cellStencil
:
    public labelList
{
    //- Owner cell (global index)
    label owner_;

    //- Location of cell centre
    vector centre_;

    //- Level of each neighbour cell
    labelList levels_;

    //- Local cells
    mutable labelList localStencil_;

    //- Local levels
    mutable labelList localLevels_;


public:

    cellStencil()
    :
        owner_(-1)
    {}

    cellStencil(const cellStencil& stencil)
    :
        labelList(stencil),
        owner_(stencil.owner_),
        centre_(stencil.centre_),
        levels_(stencil.levels_)
    {}

    cellStencil
    (
        const label own,
        const labelList& stencil,
        const labelList& levels,
        const vector& centre = vector::zero
    )
    :
        labelList(stencil),
        owner_(own),
        centre_(centre),
        levels_(levels)
    {}

    cellStencil
    (
        const label own,
        const labelList&& stencil,
        const Map<label>& levelMap,
        const vector& centre = vector::zero
    )
    :
        labelList(stencil),
        owner_(own),
        centre_(centre),
        levels_(stencil.size(), 0)
    {
        labelList& st(*this);
        forAll(st, i)
        {
            levels_[i] = levelMap[st[i]];
        }
    }

    void updateLocalStencil(const globalIndex&);

    label owner() const
    {
        return owner_;
    }

    label& owner()
    {
        return owner_;
    }

    const vector& centre() const
    {
        return centre_;
    }

    vector& centre()
    {
        return centre_;
    }

    const labelList& levels() const
    {
        return levels_;
    }

    labelList& levels()
    {
        return levels_;
    }


    const labelList& localStencil() const
    {
        return localStencil_;
    }

    labelList localStencil(const label level) const;

    void operator=(const cellStencil& stencil)
    {
        labelList::operator=(stencil);
        owner_ = stencil.owner_;
        centre_ = stencil.centre_;
        levels_ = stencil.levels_;
    }

    friend Ostream& operator<<(Ostream& os, const cellStencil&);

    friend Istream& operator>>(Istream& is, cellStencil&);
};


/*---------------------------------------------------------------------------*\
        Class extendedNLevelGlobalCellToCellStencil Declaration
\*---------------------------------------------------------------------------*/

template<class StencilType>
class extendedNLevelGlobalCellToCellStencil
:
    public MeshObject
    <
        polyMesh,
        Foam::MoveableMeshObject,
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
        mutable globalIndex gIndex_;

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
        virtual bool movePoints()
        {
            mapPtr_.clear();
            return true;
        }

        //- Update stencils. Return false if nothing changed.
        virtual bool update()
        {
            mapPtr_.clear();
            return true;
        }

        //- Collect data across all processors
        template<class Type>
        void collectData
        (
            const Field<Type>& fld,
            List<List<Type>>& stencilFld
        ) const;

        inline const globalIndex& globalNumbering() const
        {
            if (!mapPtr_.valid())
            {
                updateStencil();
            }
            return gIndex_;
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
            return stencilMap_[gIndex_.toGlobal(celli)];
        }

        virtual bool writeData(Ostream&) const
        {
            return false;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
//

#ifdef NoRepository
#include "extendedNLevelGlobalCellToCellStencil.C"
#include "extendedNLevelGlobalCellToCellStencilTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
