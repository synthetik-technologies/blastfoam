/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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
------------------------------------------------------------------------*/

#ifndef accelerationSchemeNew_H
#define accelerationSchemeNew_H

#include "AccelerationSchemeBase.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::autoPtr<Foam::accelerationScheme> Foam::accelerationScheme::New
(
    GeometricField<Type, Patch, Mesh>& field,
    autoPtr<PatchFieldSelector<Type>> selector,
    const dictionary& dict
)
{
    return AccelerationSchemeBase<Type, Patch, Mesh>::New
    (
        field,
        selector,
        dict
    );
}

#endif

// ************************************************************************* //
