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
    owner_(-1),
    maxLevel_(0),
    centre_(Zero),
    offsets_(0)
{}

Foam::cellStencil::cellStencil(const cellStencil& stencil)
:
    labelList(stencil),
    owner_(stencil.owner_),
    maxLevel_(stencil.maxLevel_),
    centre_(stencil.centre_),
    offsets_(stencil.offsets_)
{}

Foam::cellStencil::cellStencil
(
    const label own,
    const labelList& stencil,
    const labelList& levels,
    const vector& centre
)
:
    labelList(stencil.size(), -1),
    owner_(own),
    maxLevel_(0),
    centre_(centre),
    offsets_(stencil.size(), -1),
    localStencil_(),
    localOffsets_()
{
    if (stencil.size())
    {
    maxLevel_ = max(levels);
        offsets_.setSize(maxLevel_ + 2, -1);

    label ci = 0;
    for (label leveli = 0; leveli <= maxLevel_; leveli++)
    {
        offsets_[leveli] = ci;
        forAll(levels, i)
        {
        if (levels[i] == leveli)
        {
            operator[](ci++) = stencil[i];
        }
        }
    }
    }
    offsets_.last() = stencil.size();
}

Foam::cellStencil::cellStencil
(
    const label own,
    const labelList& stencil,
    const Map<label>& levelMap,
    const vector& centre
)
:
    labelList(stencil.size(), -1),
    owner_(own),
    maxLevel_(-1),
    centre_(centre),
    offsets_(stencil.size(), -1),
    localStencil_(),
    localOffsets_()
{
    if (stencil.size())
    {
    labelList levels(stencil.size(), 0);
    forAll(stencil, i)
    {
        levels[i] = levelMap[stencil[i]];
    }

    maxLevel_ = max(levels);
    offsets_.setSize(maxLevel_ + 2, -1);
    label ci = 0;
    for (label leveli = 0; leveli <= maxLevel_; leveli++)
    {
        offsets_[leveli] = ci;
        forAll(levels, i)
        {
        if (levels[i] == leveli)
        {
            operator[](ci++) = stencil[i];
        }
        }
    }
    offsets_.last() = stencil.size();
    }
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::cellStencil::updateLocalStencil
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
    localOffsets_.setSize(offsets_.size(), -1);

    if (!stencil.size())
    {
        return;
    }

    if (idx.isLocal(this->first()))
    {
        owner_ = idx.toLocal(this->first());
    }
    else
    {
        owner_ = -1;
    }

    label ni = 0;
    for (label leveli = 0; leveli <= maxLevel_; leveli++)
    {
    bool offsetSet = false;
    for (label i = offsets_[leveli]; i < offsets_[leveli+1]; i++)
    {
        if (idx.isLocal(stencil[i]))
        {
        localStencil_[ni] = idx.toLocal(stencil[i]);
            if (!offsetSet)
            {
            localOffsets_[leveli] = ni;
            offsetSet = true;
        }
        ni++;
        }
    }

        if (!offsetSet)
    {
        localOffsets_[leveli] = ni;
    }
    }

    localStencil_.resize(ni);
    localOffsets_.last() = ni;
}


Foam::SubList<Foam::label> Foam::cellStencil::operator()
(
    const label level
) const
{
    if (level > maxLevel_)
    {
    return SubList<label>(*this, 0, 0);
    }
    return SubList<label>
    (
    *this,
        offsets_[level+1] - offsets_[level],
    offsets_[level]
    );
}


Foam::SubList<Foam::label> Foam::cellStencil::operator()
(
    const label level0,
    const label level1
) const
{
    if (level0 > maxLevel_)
    {
    return SubList<label>(*this, 0, 0);
    }

    return SubList<label>
    (
    *this,
        offsets_[max(level1, maxLevel_) + 1] - offsets_[level0],
    offsets_[level0]
    );
}


Foam::SubList<Foam::label> Foam::cellStencil::localStencil
(
    const label level
) const
{
    if (!Pstream::parRun())
    {
    return operator()(level);
    }
    else if (level > maxLevel_ + 1)
    {
    return SubList<label>(localStencil_, 0, 0);
    }

    return SubList<label>
    (
    localStencil_,
        localOffsets_[level+1] - localOffsets_[level],
    localOffsets_[level]
    );
}


Foam::SubList<Foam::label> Foam::cellStencil::localStencil
(
    const label level0,
    const label level1
) const
{
    if (!Pstream::parRun())
    {
    return operator()(level0, level1);
    }
    else if (level0 > maxLevel_)
    {
    return SubList<label>(localStencil_, 0, 0);
    }
    return SubList<label>
    (
    localStencil_,
        localOffsets_[max(level1, maxLevel_) + 1] - localOffsets_[level0],
    localOffsets_[level0]
    );
}


Foam::Ostream& Foam::operator<<(Ostream& os, const cellStencil& c)
{
    const labelList& lst(c);
    os  << lst << " "
        << c.owner_ << " "
    << c.maxLevel_ << " "
        << c.offsets_ << " "
        << c.centre_ << " ";
    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, cellStencil& c)
{
    labelList& lst(c);
    is >> lst;
    is >> c.owner_;
    is >> c.maxLevel_;
    is >> c.offsets_;
    is >> c.centre_;
    return is;
}

// ************************************************************************* //
