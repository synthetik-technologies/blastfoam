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

\*---------------------------------------------------------------------------*/

#include "cellStencil.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::cellStencil::updateLocalStencil
(
    const globalIndex& idx
)
{
    const labelList& stencil(*this);
    localStencil_ = stencil;
    if (!levels_.size())
    {
        levels_.setSize(stencil.size(), 0);
    }
    localLevels_ = levels_;
    label ni = 0;
    forAll(stencil, i)
    {
        if (idx.isLocal(stencil[i]))
        {
            localStencil_[ni] = idx.toLocal(stencil[i]);
            localLevels_[ni] = levels_[i];
            ni++;
        }
    }
    localStencil_.resize(ni);
    localLevels_.resize(ni);
}


Foam::labelList Foam::cellStencil::localStencil(const label level) const
{
    labelList st(localStencil_);
    label i = 0;
    forAll(st, I)
    {
        if (localLevels_[I] == level)
        {
            st[i++] = st[I];
        }
    }
    st.resize(i);
    return st;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const cellStencil& c)
{
    const labelList& lst(c);
    os  << lst << " "
        << c.levels() << " "
        << c.owner() << " "
        << c.centre() << " ";
    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, cellStencil& c)
{
    labelList& lst(c);
    is >> lst;
    is >> c.levels();
    is >> c.owner();
    is >> c.centre();
    return is;
}

// ************************************************************************* //
