/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

inline const Foam::word& Foam::mappedPointPatchSelector::sampleRegion() const
{
    if (mappedMovingPatchPtr_)
    {
        return mappedMovingPatchPtr_->sampleRegion();
    }

    NotImplemented;
    return mappedMovingPatchPtr_->sampleRegion();
}


inline const Foam::word& Foam::mappedPointPatchSelector::samplePatch() const
{
    if (mappedMovingPatchPtr_)
    {
        return mappedMovingPatchPtr_->sampleRegion();
    }

    NotImplemented;
    return mappedMovingPatchPtr_->sampleRegion();
}


inline const Foam::word& Foam::mappedPointPatchSelector::coupleGroup() const
{
    if (mappedMovingPatchPtr_)
    {
        return mappedMovingPatchPtr_->coupleGroup();
    }

    NotImplemented;
    return mappedMovingPatchPtr_->coupleGroup();
}


inline Foam::label Foam::mappedPointPatchSelector::sampleSize() const
{
    if (mappedMovingPatchPtr_)
    {
        return mappedMovingPatchPtr_->sampleSize();
    }

    NotImplemented;
    return mappedMovingPatchPtr_->sampleSize();
}


inline bool Foam::mappedPointPatchSelector::sameRegion() const
{
    if (mappedMovingPatchPtr_)
    {
        return mappedMovingPatchPtr_->sameRegion();
    }

    NotImplemented;
    return mappedMovingPatchPtr_->sameRegion();
}


inline const Foam::mapDistribute& Foam::mappedPointPatchSelector::map() const
{
    if (mappedMovingPatchPtr_)
    {
        return mappedMovingPatchPtr_->map();
    }

    NotImplemented;
    return mappedMovingPatchPtr_->map();
}


// ************************************************************************* //
