/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------

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

inline const Foam::word& Foam::mappedPatchSelector::sampleRegion() const
{
    if (mappedPatchPtr_)
    {
        return mappedPatchPtr_->sampleRegion();
    }

    return mappedMovingPatchPtr_->sampleRegion();
}


inline const Foam::word& Foam::mappedPatchSelector::samplePatch() const
{
    if (mappedPatchPtr_)
    {
        return mappedPatchPtr_->sampleRegion();
    }

    return mappedMovingPatchPtr_->sampleRegion();
}


inline const Foam::word& Foam::mappedPatchSelector::coupleGroup() const
{
    if (mappedPatchPtr_)
    {
        return mappedPatchPtr_->coupleGroup();
    }

    return mappedMovingPatchPtr_->coupleGroup();
}


inline Foam::label Foam::mappedPatchSelector::sampleSize() const
{
    if (mappedPatchPtr_)
    {
        return mappedPatchPtr_->sampleSize();
    }

    return mappedMovingPatchPtr_->sampleSize();
}


inline const Foam::mapDistribute& Foam::mappedPatchSelector::map() const
{
    if (mappedPatchPtr_)
    {
        return mappedPatchPtr_->map();
    }

    return mappedMovingPatchPtr_->map();
}


// ************************************************************************* //
