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

\*---------------------------------------------------------------------------*/

#include "fluxScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> fluxScheme::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& f,
    const word& name
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;
    label nCmpts = pTraits<Type>::nComponents;

    fieldType fOwn(fvc::interpolate(f, own_(), scheme(name)));
    fieldType fNei(fvc::interpolate(f, nei_(), scheme(name)));

    tmp<fieldType> tmpf
    (
        new fieldType
        (
            IOobject
            (
                f.name() + "f",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensioned<Type>("0", f.dimensions(), pTraits<Type>::zero)
        )
    );
    fieldType& ff = tmpf.ref();

    forAll(fOwn, facei)
    {
        for (label i = 0; i < nCmpts; i++)
        {
            setComponent(ff[facei], i) = interpolate
            (
                component(fOwn[facei], i),
                component(fNei[facei], i),
                facei
            );
        }
    }

    forAll(f.boundaryField(), patchi)
    {
        forAll(f.boundaryField()[patchi], facei)
        {
            for (label i = 0; i < nCmpts; i++)
            {
                setComponent(ff.boundaryFieldRef()[patchi][facei], i) =
                    interpolate
                    (
                        component(fOwn.boundaryField()[patchi][facei], i),
                        component(fNei.boundaryField()[patchi][facei], i),
                        facei, patchi
                    );
            }
        }
    }
    return tmpf;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
