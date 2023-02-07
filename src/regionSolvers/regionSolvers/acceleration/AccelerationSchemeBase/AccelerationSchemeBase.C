/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
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

#include "AccelerationSchemeBase.H"
#include "UautoPtr.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::AccelerationSchemeBase<Type, Patch, Mesh>::AccelerationSchemeBase
(
    const word& type,
    GeometricField<Type, Patch, Mesh>& field,
    const label patchi,
    const dictionary& dict
)
:
    accelerationScheme(type, patchi, dict),
    field_(field)
{
    // Make sure boundaries are actually fixed
    if (!field_.boundaryField()[patchi_].fixesValue())
    {
        FatalErrorInFunction
            << "Trying to relax a non-fixed boundary for "
            << field_.name() << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::AccelerationSchemeBase<Type, Patch, Mesh>::~AccelerationSchemeBase()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
const Foam::dictionary&
Foam::AccelerationSchemeBase<Type, Patch, Mesh>::coeffDict
(
    const dictionary& dict
) const
{
    UautoPtr<const dictionary> regionDict(dict.subDictPtr(field_.mesh().name()));
    UautoPtr<const dictionary> fieldDict;
    UautoPtr<const dictionary> patchDict;

    if (regionDict.valid())
    {
        fieldDict.set(regionDict->subDictPtr(field_.name()));
    }
    if (fieldDict.valid())
    {
        patchDict.set
        (
            fieldDict->subDictPtr(field_.mesh().boundary()[patchi_].name())
        );
    }

    if (patchDict.valid() && patchDict->found(accelerationScheme::typeName))
    {
        return patchDict->optionalSubDict(type_ + "Coeffs");
    }
    else if (fieldDict.valid()&& fieldDict->found(accelerationScheme::typeName))
    {
        return fieldDict->optionalSubDict(type_ + "Coeffs");
    }
    else if
    (
        regionDict.valid()
     && regionDict->found(accelerationScheme::typeName)
    )
    {
        return regionDict->optionalSubDict(type_ + "Coeffs");
    }
    return dict.optionalSubDict(type_ + "Coeffs");
}


template<class Type, template<class> class Patch, class Mesh>
void Foam::AccelerationSchemeBase<Type, Patch, Mesh>::updateError()
{
    error_ = Zero;
    const GeometricField<Type, Patch, Mesh>& fieldPrev = field_.prevIter();

    const Field<Type>& pfield =
        dynamicCast<const Field<Type>>(field_.boundaryField()[patchi_]);
    const Field<Type>& pfieldPrev =
        dynamicCast<const Field<Type>>(fieldPrev.boundaryField()[patchi_]);
    forAll(pfield, i)
    {
        error_ += magSqr(pfield[i] - pfieldPrev[i]);
    }

    reduce(error_, sumOp<scalar>());
    error_ = sqrt(error_);

    if (initialError_ < 0)
    {
        initialError_ = error_;
    }
}


// ************************************************************************* //
