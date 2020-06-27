/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "ensightPTraits.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::ensightPTraits<Foam::scalar>::typeName =
    Foam::pTraits<Foam::scalar>::typeName;

const char* const Foam::ensightPTraits<Foam::vector>::typeName =
    Foam::pTraits<Foam::vector>::typeName;

const char* const Foam::ensightPTraits<Foam::sphericalTensor>::typeName =
    Foam::pTraits<Foam::scalar>::typeName;

const char* const Foam::ensightPTraits<Foam::symmTensor>::typeName =
    "tensor symm";

const char* const Foam::ensightPTraits<Foam::tensor>::typeName =
    "tensor asym";


// ************************************************************************* //
