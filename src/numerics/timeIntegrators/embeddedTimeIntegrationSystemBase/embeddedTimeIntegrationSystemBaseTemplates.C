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

#include "embeddedTimeIntegrationSystemBase.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<template<class> class ListType, class Type>
Type Foam::embeddedTimeIntegrationSystemBase::calcError
(
    const ListType<Type>& fList
) const
{
    // List of error coefficients
    const List<scalar>& errorCoeffs = timeInt_->bs().last();

    // Remove old steps
    Type error(Zero);
    forAll(errorCoeffs, stepi)
    {
        label fi = timeInt_->getDeltaIndex(stepi);
        if (fi != -1 && errorCoeffs[fi] != 0)
        {
            error += errorCoeffs[stepi]*fList[fi];
        }
    }
    return error;
}

template<class Type>
Foam::scalar Foam::embeddedTimeIntegrationSystemBase::normaliseError
(
    const Type& y0,
    const Type& y,
    const Type& err,
    const scalar relTol,
    const scalar absTol
) const
{
    // Calculate the maximum error
    Type tol
    (
        pTraits<Type>::one*absTol
      + relTol*max(cmptMag(y0), cmptMag(y))
    );
    return cmptMax(cmptDivide(cmptMag(err), tol));
}
// ************************************************************************* //
