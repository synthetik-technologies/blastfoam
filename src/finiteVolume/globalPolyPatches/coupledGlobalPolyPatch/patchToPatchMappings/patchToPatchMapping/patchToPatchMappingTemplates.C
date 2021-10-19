/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "patchToPatchMapping.H"

namespace Foam
{

// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //

template<class Type>
tmp<Field<Type>>  patchToPatchMapping::transferFaces
(
    const standAlonePatch& fromPatch,
    const Field<Type>& fromField
) const
{
    if (&fromPatch == &zoneA())
    {
        tmp<Field<Type>> toFieldTmp(new Field<Type>(zoneB().size()));
        transferFaces
        (
            fromPatch,
            zoneB(),
            fromField,
            toFieldTmp.ref()
        );
        return toFieldTmp;
    }
    else if (&fromPatch == &zoneB())
    {
        tmp<Field<Type>> toFieldTmp(new Field<Type>(zoneA().size()));
        transferFaces
        (
            fromPatch,
            zoneA(),
            fromField,
            toFieldTmp.ref()
        );
        return toFieldTmp;
    }
    else
    {
        FatalErrorInFunction
            << "Patch does not belong to this mapping"
            << abort(FatalError);
    }
    return tmp<Field<Type>>();
}


template<class Type>
tmp<Field<Type>>  patchToPatchMapping::transferFaces
(
    const standAlonePatch& fromPatch,
    const tmp<Field<Type>>& fromField
) const
{
    return transferFaces(fromPatch, fromField());
}


template<class Type>
tmp<Field<Type>>  patchToPatchMapping::transferPoints
(
    const standAlonePatch& fromPatch,
    const Field<Type>& fromField
) const
{
    if (&fromPatch == &zoneA())
    {
        tmp<Field<Type>> toFieldTmp(new Field<Type>(zoneB().nPoints()));
        transferPoints
        (
            fromPatch,
            zoneB(),
            fromField,
            toFieldTmp.ref()
        );
        return toFieldTmp;
    }
    else if (&fromPatch == &zoneB())
    {
        tmp<Field<Type>> toFieldTmp(new Field<Type>(zoneA().nPoints()));
        transferPoints
        (
            fromPatch,
            zoneA(),
            fromField,
            toFieldTmp.ref()
        );
        return toFieldTmp;
    }
    else
    {
        FatalErrorInFunction
            << "Patch does not belong to this mapping"
            << abort(FatalError);
    }
    return tmp<Field<Type>>();
}

template<class Type>
tmp<Field<Type>>  patchToPatchMapping::transferPoints
(
    const standAlonePatch& fromPatch,
    const tmp<Field<Type>>& fromField
) const
{
    return transferPoints(fromPatch, fromField());
}

// ************************************************************************* //


} // End namespace Foam

// ************************************************************************* //
