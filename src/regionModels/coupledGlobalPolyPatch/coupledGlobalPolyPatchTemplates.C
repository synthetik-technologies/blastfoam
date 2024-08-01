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

#include "coupledGlobalPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledGlobalPolyPatch::pointInterpolate
(
    const Field<Type>& pField,
    const Field<Type>& unmapped
) const
{
    if (pField.size() != this->physicalPatch().nPoints())
    {
        FatalErrorInFunction
            << "Patch field does not correspond to patch points." << nl
            << "Patch size: " << this->physicalPatch().nPoints() << nl
            << "Field size: " << pField.size()
            << abort(FatalError);
        return tmp<Field<Type>>();
    }

    // TODO handle point interpolation
    const patchToPatchMapping& mapper = patchToPatchInterpolator(true);
    if (isSrc_)
    {
        return mapper.transferPointsToTgt(pField, unmapped);
    }
    else
    {
        return mapper.transferPointsToSrc(pField, unmapped);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledGlobalPolyPatch::pointInterpolate
(
    const tmp<Field<Type>>& tpField,
    const Field<Type>& unmapped
) const
{
    return pointInterpolate(tpField(), unmapped);
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::coupledGlobalPolyPatch::faceInterpolate
(
    const Field<Type>& fField,
    const Field<Type>& unmapped
) const
{
    if (fField.size() != this->physicalPatch().size())
    {
        FatalErrorInFunction
            << "Patch field does not correspond to patch faces." << nl
            << "Patch size: " << this->physicalPatch().size() << nl
            << "Field size: " << fField.size()
            << abort(FatalError);
        return tmp<Field<Type>>();
    }

    const patchToPatchMapping& mapper = patchToPatchInterpolator(false);
    if (isSrc_)
    {
        return
            mapper.transferFacesToTgt
            (
                fField,
                unmapped
            );
    }
    else
    {
        return
            mapper.transferFacesToSrc
            (
                fField,
                unmapped
            );
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::coupledGlobalPolyPatch::faceInterpolate
(
    const tmp<Field<Type>>& fField,
    const Field<Type>& unmapped
) const
{

    return faceInterpolate(fField(), unmapped);
}


template<class Type>
void Foam::coupledGlobalPolyPatch::setUnmappedFace
(
    Field<Type>& ff,
    const Type& unmapped
) const
{
    UIndirectList<Type>(ff, unmappedFaces_) = unmapped;
}


template<class Type>
void Foam::coupledGlobalPolyPatch::setUnmappedFace
(
    Field<Type>& ff,
    const Field<Type>& unmapped
) const
{
    UIndirectList<Type>(ff, unmappedFaces_) =
        UIndirectList<Type>(unmapped, unmappedFaces_);
}


template<class Type>
void Foam::coupledGlobalPolyPatch::setUnmappedFace
(
    Field<Type>& ff,
    const tmp<Field<Type>>& tunmapped
) const
{
    setUnmappedFace(ff, tunmapped());
}


template<class Type>
void Foam::coupledGlobalPolyPatch::setUnmappedPoint
(
    Field<Type>& pf,
    const Type& unmapped
) const
{
    UIndirectList<Type>(pf, unmappedPoints_) = unmapped;
}


template<class Type>
void Foam::coupledGlobalPolyPatch::setUnmappedPoint
(
    Field<Type>& pf,
    const Field<Type>& unmapped
) const
{
    UIndirectList<Type>(pf, unmappedPoints_) =
        UIndirectList<Type>(unmapped, unmappedPoints_);
}


template<class Type>
void Foam::coupledGlobalPolyPatch::setUnmappedPoint
(
    Field<Type>& pf,
    const tmp<Field<Type>>& tunmapped
) const
{
    setUnmappedFace(pf, tunmapped());
}


// ************************************************************************* //
