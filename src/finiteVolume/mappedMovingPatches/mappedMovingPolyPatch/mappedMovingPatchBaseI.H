/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-08-21 Synthetik Applied Technology: Mapping of patches
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


inline Foam::tmp<Foam::vectorField> Foam::mappedMovingPatchBase::offsets() const
{
    const polyMesh& sMesh(sampleMesh());
    const polyPatch& sPatch(samplePolyPatch());
    if (sMesh.foundObject<pointVectorField>(displacementFieldName_))
    {
        return
            sMesh.lookupObject<pointVectorField>
            (
                displacementFieldName_
            ).boundaryField()[sPatch.index()].patchInternalField();
    }

    if (debug)
    {
        Info<< "Offsets are not set for " << samplePatch() << endl;
    }

    return tmp<vectorField>(new vectorField(sPatch.localPoints().size(), Zero));
}


inline const Foam::word& Foam::mappedMovingPatchBase::displacementFieldName() const
{
    return displacementFieldName_;
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
