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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellStencil::cellStencil()
:
    labelList(0),
    localOwner_(-1),
    centre_(Zero),
    localStencil_()
{}

Foam::cellStencil::cellStencil(const cellStencil& stencil)
:
    labelList(stencil),
    localOwner_(stencil.localOwner_),
    centre_(stencil.centre_),
    localStencil_()
{}

Foam::cellStencil::cellStencil
(
    const label own,
    const labelList& stencil,
    const vector& centre
)
:
    labelList(stencil),
    localOwner_(own),
    centre_(centre),
    localStencil_()
{
    if (stencil.size())
    {
        if (stencil[0] != own)
        {
            label i = 0;
            this->setSize(stencil.size() + 1);
            operator[](i++) = own;
            forAll(stencil, si)
            {
                operator[](i++) = stencil[si];
            }
        }
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::cellStencil::update
(
    const globalIndex& idx
)
{
    if (!Pstream::parRun())
    {
        return;
    }

    const labelList& stencil(*this);
    localStencil_.setSize(stencil.size(), -1);

    if (!this->size())
    {
        localOwner_ = -1;
        return;
    }

    if (idx.isLocal(this->first()))
    {
        localOwner_ = idx.toLocal(this->first());
    }
    else
    {
        localOwner_ = -1;
    }

    if (!stencil.size())
    {
        return;
    }

    label li = 0;
    forAll(*this, i)
    {
        if (idx.isLocal(labelList::operator[](i)))
        {
            localStencil_[li++] = idx.toLocal(labelList::operator[](i));
        }
    }
    localStencil_.resize(li);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const cellStencil& c)
{
    const labelList& lst(c);
    os  << lst << " "
        << c.localOwner_ << " "
        << c.centre_ << " ";
    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, cellStencil& c)
{
    labelList& lst(c);
    is >> lst;
    is >> c.localOwner_;
    is >> c.centre_;
    return is;
}

// ************************************************************************* //
