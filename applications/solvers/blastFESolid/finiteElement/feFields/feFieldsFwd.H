/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

InClass
    Foam::feFields

Description

SourceFiles
    fetFields.C

\*---------------------------------------------------------------------------*/

#ifndef feFieldsFwd_H
#define feFieldsFwd_H

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class feMesh;

template<class Type>
class fePatchField;

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricField;

template<class Type>
using FEField = GeometricField<Type, fePatchField, feMesh>;

typedef GeometricField<scalar, fePatchField, feMesh> feScalarField;
typedef GeometricField<vector, fePatchField, feMesh> feVectorField;
typedef GeometricField<sphericalTensor, fePatchField, feMesh>
    feSphericalTensorField;
typedef GeometricField<symmTensor, fePatchField, feMesh>
    feSymmTensorField;
typedef GeometricField<tensor, fePatchField, feMesh> feTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
