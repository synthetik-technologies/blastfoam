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
    Foam::cellStencil

Description
    Basic class for a cell stencil

SourceFiles
    cellStencil.C

\*---------------------------------------------------------------------------*/

#ifndef cellStencil_H
#define cellStencil_H

#include "labelList.H"
#include "vector.H"
#include "Map.H"
#include "globalIndex.H"

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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
