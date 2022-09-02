/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
03-12-2021 Synthetik Applied Technologies : Added Function3
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

\*---------------------------------------------------------------------------*/

#include "Function3Evaluate.H"

// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::evaluate
(
    GeometricField<Type, PatchField, GeoMesh>& result,
    const Function3<Type>& func,
    const GeometricField<Type, PatchField, GeoMesh>& x,
    const GeometricField<Type, PatchField, GeoMesh>& y,
    const GeometricField<Type, PatchField, GeoMesh>& z
)
{
    result.primitiveFieldRef() = func.value(x(), y(), z());

    typename GeometricField<Type, PatchField, GeoMesh>::Boundary& bresult =
        result.boundaryFieldRef();

    const typename GeometricField<Type, PatchField, GeoMesh>::Boundary& bx =
        x.boundaryField();

    const typename GeometricField<Type, PatchField, GeoMesh>::Boundary& by =
        y.boundaryField();

    const typename GeometricField<Type, PatchField, GeoMesh>::Boundary& bz =
        z.boundaryField();

    forAll(bresult, patchi)
    {
        bresult[patchi] = func.value
        (
            bx[patchi],
            by[patchi],
            bz[patchi]
        );
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>> Foam::evaluate
(
    const Function3<Type>& func,
    const dimensionSet& dims,
    const GeometricField<Type, PatchField, GeoMesh>& x,
    const GeometricField<Type, PatchField, GeoMesh>& y,
    const GeometricField<Type, PatchField, GeoMesh>& z
)
{
    tmp<GeometricField<Type, PatchField, GeoMesh>> tresult
    (
        GeometricField<Type, PatchField, GeoMesh>::New
        (
            func.name()
          + '(' + x.name() + ',' + y.name() + ',' + z.name() ')',
            x.mesh(),
            dims
        )
    );

    evaluate(tresult.ref(), func, x, y, z);

    return tresult;
}


// ************************************************************************* //
