/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
01-11-2022 Synthetik Applied Technologies : Added Function4
-------------------------------------------------------------------------------
License
    This file is a derived work of OpenFOAM.

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

Class
    Foam::Function4Evaluate

Description
    Global functions to evaluate Function4 of GeometricFields

SourceFiles
    Function4Evaluate.C

\*---------------------------------------------------------------------------*/

#ifndef Function4Evaluate_H
#define Function4Evaluate_H

#include "Function4.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void evaluate
(
    GeometricField<Type, PatchField, GeoMesh>& result,
    const Function4<Type>& func,
    const GeometricField<Type, PatchField, GeoMesh>& x,
    const GeometricField<Type, PatchField, GeoMesh>& y,
    const GeometricField<Type, PatchField, GeoMesh>& z
);

template<class Type, template<class> class PatchField, class GeoMesh>
void evaluate
(
    GeometricField<Type, PatchField, GeoMesh>& result,
    const Function4<Type>& func,
    const scalar t,
    const GeometricField<Type, PatchField, GeoMesh>& x,
    const GeometricField<Type, PatchField, GeoMesh>& y,
    const GeometricField<Type, PatchField, GeoMesh>& z
);

template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> evaluate
(
    const Function4<Type>& func,
    const dimensionSet& dims,
    const GeometricField<Type, PatchField, GeoMesh>& x,
    const GeometricField<Type, PatchField, GeoMesh>& y,
    const GeometricField<Type, PatchField, GeoMesh>& z
);


template<class Type, template<class> class PatchField, class GeoMesh>
tmp<GeometricField<Type, PatchField, GeoMesh>> evaluate
(
    const Function4<Type>& func,
    const dimensionSet& dims,
    const scalar t,
    const GeometricField<Type, PatchField, GeoMesh>& x,
    const GeometricField<Type, PatchField, GeoMesh>& y,
    const GeometricField<Type, PatchField, GeoMesh>& z
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Function4Evaluate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
