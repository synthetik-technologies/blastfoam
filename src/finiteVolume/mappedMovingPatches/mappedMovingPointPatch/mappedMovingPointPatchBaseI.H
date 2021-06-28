/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-08-21 Synthetik Applied Technology: Mapping of point patches
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


inline const Foam::word& Foam::mappedMovingPointPatchBase::sampleRegion() const
{
    return mpp_.sampleRegion();
}


inline const Foam::word& Foam::mappedMovingPointPatchBase::samplePatch() const
{
    return mpp_.samplePatch();
}


inline const Foam::word& Foam::mappedMovingPointPatchBase::coupleGroup() const
{
    return mpp_.coupleGroup();
}


inline Foam::label Foam::mappedMovingPointPatchBase::sampleSize() const
{
    return samplePolyPatch().nPoints();
}


inline Foam::tmp<Foam::vectorField> Foam::mappedMovingPointPatchBase::offsets() const
{
    const polyMesh& sMesh(sampleMesh());
    const polyPatch sPatch(samplePolyPatch());

    if (sMesh.foundObject<pointVectorField>(displacementFieldName_))
    {
        tmp<vectorField> disp
        (
            sMesh.lookupObject<pointVectorField>
            (
                displacementFieldName_
            ).boundaryField()[sPatch.index()].patchInternalField()
        );
        return disp;
    }

    if (debug)
    {
        Info<< "Offsets are not set for " << sPatch.name() << endl;
    }

    return tmp<vectorField>(new vectorField(sPatch.size(), Zero));
}


inline const Foam::mapDistribute& Foam::mappedMovingPointPatchBase::map() const
{
    if (mapPtr_.empty())
    {
        calcMapping();
    }

    return mapPtr_();
}

// ************************************************************************* //
