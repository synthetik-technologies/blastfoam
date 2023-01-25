/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

SourceFiles
    fvcInterpolate.C

\*---------------------------------------------------------------------------*/


#ifndef fvcInterpolate_H
#define fvcInterpolate_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "pointFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace pointFieldOps
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void pushUntransformedData(List<Type>& pointData, const polyMesh& mesh);

template<class Type>
void addSeparated(GeometricField<Type, pointPatchField, pointMesh>& pf) ;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pointFieldOps

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Volume to point interpolation with gradient reconstruction
template<class Type>
void volPointInterpolate
(
    const GeometricField<Type, fvPatchField, volMesh>&,
    const GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >&,
    GeometricField<Type, pointPatchField, pointMesh>&,
    const bool interpolateBoundary = false
);

//- Surface to volume interpolation
template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> surfVolInterpolate
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>&
);


//- Point to surface interpolation
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> pointSurfInterpolate
(
    const GeometricField<Type, pointPatchField, pointMesh>&
);


//- Point to volume interpolation
template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>> pointVolInterpolate
(
    const GeometricField<Type, pointPatchField, pointMesh>&
);

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fvcInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
