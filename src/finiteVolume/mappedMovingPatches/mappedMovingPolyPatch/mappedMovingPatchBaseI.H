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

inline const Foam::word& Foam::mappedMovingPatchBase::sampleRegion() const
{
    if (sampleRegion_.empty())
    {
        if (!coupleGroup_.valid())
        {
            FatalErrorInFunction
                << "Supply either a regionName or a coupleGroup"
                << " for patch " << patch_.name()
                << " in region " << patch_.boundaryMesh().mesh().name()
                << exit(FatalError);
        }

        // Try and use patchGroup to find samplePatch and sampleRegion
        label samplePatchID = coupleGroup_.findOtherPatchID
        (
            patch_,
            sampleRegion_
        );

        samplePatch_ = sampleMesh().boundaryMesh()[samplePatchID].name();
    }
    return sampleRegion_;
}


inline const Foam::word& Foam::mappedMovingPatchBase::samplePatch() const
{
    if (samplePatch_.empty())
    {
        if (!coupleGroup_.valid())
        {
            FatalErrorInFunction
                << "Supply either a patchName or a coupleGroup"
                << " for patch " << patch_.name()
                << " in region " << patch_.boundaryMesh().mesh().name()
                << exit(FatalError);
        }

        // Try and use patchGroup to find samplePatch and sampleRegion
        label samplePatchID = coupleGroup_.findOtherPatchID
        (
            patch_,
            sampleRegion_
        );

        samplePatch_ = sampleMesh().boundaryMesh()[samplePatchID].name();
    }
    return samplePatch_;
}


inline const Foam::word& Foam::mappedMovingPatchBase::coupleGroup() const
{
    return coupleGroup_.name();
}


inline Foam::label Foam::mappedMovingPatchBase::sampleSize() const
{
    return samplePolyPatch().size();
}


inline void Foam::mappedMovingPatchBase::setOffsets
(
    const volVectorField& D
)
{
    if (DPtr_)
    {
        return;
    }

    if (debug)
    {
        Info<< "Setting offsets for " << patch_.name() << endl;
    }
    DPtr_ = &D;
}


inline const Foam::mapDistribute& Foam::mappedMovingPatchBase::map() const
{
    if (mapPtr_.empty())
    {
        calcMapping();
    }

    return mapPtr_();
}


// ************************************************************************* //
