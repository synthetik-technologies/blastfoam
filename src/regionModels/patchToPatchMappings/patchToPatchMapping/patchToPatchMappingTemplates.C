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
#include "patchToPatchTools.H"

// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::patchToPatchMapping::transferFacesToTgt
(
    const Field<Type>& srcField,
    const Field<Type>& unmapped
) const
{
    if (!isNull(unmapped))
    {
        return patchToPatchTools::interpolate
        (
            localSrcFacesToTgt_,
            tgtFaceWeights_,
            srcFacesMapPtr_,
            srcField,
            unmapped
        );
    }
    else
    {
        return patchToPatchTools::interpolate
        (
            localSrcFacesToTgt_,
            tgtFaceWeights_,
            srcFacesMapPtr_,
            srcField
        );
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::patchToPatchMapping::transferFacesToTgt
(
    const tmp<Field<Type>>& tsrcField,
    const Field<Type>& unmapped
) const
{
    return transferFacesToTgt(tsrcField(), unmapped);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>  Foam::patchToPatchMapping::transferFacesToSrc
(
    const Field<Type>& tgtField,
    const Field<Type>& unmapped
) const
{
    if (!isNull(unmapped))
    {
        return patchToPatchTools::interpolate
        (
            localTgtFacesToSrc_,
            srcFaceWeights_,
            tgtFacesMapPtr_,
            tgtField,
            unmapped
        );
    }
    else
    {
        return patchToPatchTools::interpolate
        (
            localTgtFacesToSrc_,
            srcFaceWeights_,
            tgtFacesMapPtr_,
            tgtField
        );
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>  Foam::patchToPatchMapping::transferFacesToSrc
(
    const tmp<Field<Type>>& ttgtField,
    const Field<Type>& unmapped
) const
{
    return transferFacesToSrc(ttgtField(), unmapped);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::patchToPatchMapping::transferPointsToTgt
(
    const Field<Type>& srcField,
    const Field<Type>& unmapped
) const
{
    if (!isNull(unmapped))
    {
        return patchToPatchTools::interpolate
        (
            localSrcPointsToTgt_,
            tgtPointWeights_,
            srcPointsMapPtr_,
            srcField,
            unmapped
        );
    }
    else
    {
        return patchToPatchTools::interpolate
        (
            localSrcPointsToTgt_,
            tgtPointWeights_,
            srcPointsMapPtr_,
            srcField
        );
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::patchToPatchMapping::transferPointsToTgt
(
    const tmp<Field<Type>>& tsrcField,
    const Field<Type>& unmapped
) const
{
    return transferPointsToTgt(tsrcField(), unmapped);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>  Foam::patchToPatchMapping::transferPointsToSrc
(
    const Field<Type>& tgtField,
    const Field<Type>& unmapped
) const
{
    if (!isNull(unmapped))
    {
        return patchToPatchTools::interpolate
        (
            localTgtPointsToSrc_,
            srcPointWeights_,
            tgtPointsMapPtr_,
            tgtField,
            unmapped
        );
    }
    else
    {
        return patchToPatchTools::interpolate
        (
            localTgtPointsToSrc_,
            srcPointWeights_,
            tgtPointsMapPtr_,
            tgtField
        );
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>  Foam::patchToPatchMapping::transferPointsToSrc
(
    const tmp<Field<Type>>& ttgtField,
    const Field<Type>& unmapped
) const
{
    return transferPointsToSrc(ttgtField(), unmapped);
}


// ************************************************************************* //
