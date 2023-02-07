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

#include "accelerationScheme.H"
#include "UautoPtr.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(accelerationScheme, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::accelerationScheme::accelerationScheme
(
    const word& type,
    const label patchi,
    const dictionary& dict
)
:
    regionName_("unknown"),
    fieldName_("unknown"),
    patchName_("unknown"),
    type_(type),
    patchi_(patchi),
    initialError_(-1.0),
    error_(great)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::accelerationScheme::~accelerationScheme()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::accelerationScheme::schemesDict
(
    const dictionary& dict,
    const word& regionName,
    const word& fieldName,
    const word& patchName
)
{
    UautoPtr<const dictionary> regionDict(dict.subDictPtr(regionName));
    UautoPtr<const dictionary> fieldDict;
    UautoPtr<const dictionary> patchDict;

    if (regionDict.valid())
    {
        fieldDict.set(regionDict->subDictPtr(fieldName));
    }
    if (fieldDict.valid())
    {
        patchDict.set
        (
            fieldDict->subDictPtr(patchName)
        );
    }

    if (patchDict.valid() && patchDict->found(accelerationScheme::typeName))
    {
        return patchDict();
    }
    else if (fieldDict.valid()&& fieldDict->found(accelerationScheme::typeName))
    {
        return fieldDict();
    }
    else if
    (
        regionDict.valid()
     && regionDict->found(accelerationScheme::typeName)
    )
    {
        return regionDict();
    }
    return dict;
}


const Foam::dictionary& Foam::accelerationScheme::coeffDict
(
    const dictionary& dict
) const
{
    return
        schemesDict
        (
            dict,
            regionName_,
            fieldName_,
            patchName_
        ).optionalSubDict(type_ + "Coeffs");
}

template<>
Foam::scalar Foam::accelerationScheme::sumDotDot
(
    const UList<scalar>& f
) const
{
    scalar sumf = Zero;
    forAll(f, i)
    {
        sumf += sqr(f[i]);
    }
    return returnReduce(sumf, sumOp<scalar>());
}


template<>
Foam::scalar Foam::accelerationScheme::sumDotDot
(
    const UList<scalar>& f1,
    const UList<scalar>& f2
) const
{
    scalar sumf12 = Zero;
    forAll(f1, i)
    {
        sumf12 += f1[i]*f2[i];
    }
    return returnReduce(sumf12, sumOp<scalar>());
}


template<>
Foam::scalar Foam::accelerationScheme::dotDot
(
    const scalar& f
) const
{
    return f*f;
}


template<>
Foam::scalar Foam::accelerationScheme::dotDot
(
    const scalar& f1,
    const scalar& f2
) const
{
    return f1*f2;
}


void Foam::accelerationScheme::clear()
{
    initialError_ = -1;
    error_ = great;
}


// ************************************************************************* //
