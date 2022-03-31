/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

InClass
    Foam::mechanicalLaw

\*---------------------------------------------------------------------------*/

#include "mechanicalLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, template<class> class PatchField, class GeoMesh>
autoPtr<GeometricField<Type, PatchField, GeoMesh>> mechanicalLaw::makeField
(
    const word& name,
    const dimensioned<Type>& dt,
    const IOobject::readOption rOpt,
    const IOobject::writeOption wOpt
) const
{
    return autoPtr<GeometricField<Type, PatchField, GeoMesh>>
    (
        new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                name,
                mesh_.time().timeName(mesh_.time().startTime().value()),
                mesh_,
                rOpt,
                wOpt
            ),
            mesh_,
            dt
        )
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
autoPtr<GeometricField<Type, PatchField, GeoMesh>>
mechanicalLaw::makeTypeField
(
    const word& name,
    const dimensioned<Type>& dt,
    const IOobject::readOption rOpt,
    const IOobject::writeOption wOpt
) const
{
    return autoPtr<GeometricField<Type, PatchField, GeoMesh>>
    (
        new GeometricField<Type, PatchField, GeoMesh>
        (
            IOobject
            (
                IOobject::groupName(name, name_),
                mesh_.time().timeName(mesh_.time().startTime().value()),
                mesh_,
                rOpt,
                wOpt
            ),
            mesh_,
            dt
        )
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
autoPtr<GeometricField<Type, PatchField, GeoMesh>>
mechanicalLaw::makeField
(
    const IOobject& io,
    const dimensioned<Type>& dt
) const
{
    return autoPtr<GeometricField<Type, PatchField, GeoMesh>>
    (
        new GeometricField<Type, PatchField, GeoMesh>
        (
            io,
            mesh_,
            dt
        )
    );
}


template<class FieldType, class CoeffType>
void mechanicalLaw::replacePlaneStress
(
    FieldType& epsilon,
    const CoeffType& coeff,
    const FieldType& sigma,
    const FieldType& epsilonP
) const
{
    if (usePlaneStress_)
    {
        epsilon.replace
        (
            planeStressDir_,
           -coeff
           *(
               sigma.component(nonPlaneStressDirs_[0])
             + sigma.component(nonPlaneStressDirs_[1])
            )
          - (
                epsilonP.component(nonPlaneStressDirs_[0])
              + epsilonP.component(nonPlaneStressDirs_[1])
            )
        );
    }
}


template<class FieldType, class CoeffType>
void mechanicalLaw::replacePlaneStress
(
    FieldType& epsilon,
    const CoeffType& coeff,
    const FieldType& sigma
) const
{
    if (usePlaneStress_)
    {
        epsilon.replace
        (
            planeStressDir_,
           -coeff
           *(
               sigma.component(nonPlaneStressDirs_[0])
             + sigma.component(nonPlaneStressDirs_[1])
            )
        );
    }
}


template<template<class> class PatchField, class Mesh, class CoeffType>
void mechanicalLaw::updateEpsilon
(
    GeometricField<symmTensor, PatchField, Mesh>& epsilon,
    const CoeffType& coeff,
    const GeometricField<symmTensor, PatchField, Mesh>& sigma,
    const GeometricField<symmTensor, PatchField, Mesh>& epsilonP
) const
{
    word ext = word::null;
    if (epsilon.size() == mesh_.nInternalFaces())
    {
        ext = "f";
    }
    // Calculate total strain
    if (incremental())
    {
        // Lookup gradient of displacement increment
        const GeometricField<tensor, PatchField, Mesh>& gradDD =
            mesh().lookupObject
            <
                GeometricField<tensor, PatchField, Mesh>
            >("grad(DD)" + ext);

        epsilon = epsilon.oldTime() + symm(gradDD);
    }
    else
    {
        // Lookup gradient of displacement
        const GeometricField<tensor, PatchField, Mesh>& gradD =
            mesh().lookupObject
            <
                GeometricField<tensor, PatchField, Mesh>
            >("grad(D)" + ext);
        epsilon = symm(gradD);
    }
    replacePlaneStress(epsilon, coeff, sigma, epsilonP);
}


template<template<class> class PatchField, class Mesh, class CoeffType>
void mechanicalLaw::updateEpsilon
(
    GeometricField<symmTensor, PatchField, Mesh>& epsilon,
    const CoeffType& coeff,
    const GeometricField<symmTensor, PatchField, Mesh>& sigma
) const
{
    word ext = word::null;
    if (epsilon.size() == mesh_.nInternalFaces())
    {
        ext = "f";
    }

    // Calculate total strain
    if (incremental())
    {
        // Lookup gradient of displacement increment
        const GeometricField<tensor, PatchField, Mesh>& gradDD =
            mesh().lookupObject
            <
                GeometricField<tensor, PatchField, Mesh>
            >("grad(DD)" + ext);

        epsilon = epsilon.oldTime() + symm(gradDD);
    }
    else
    {
        // Lookup gradient of displacement
        const GeometricField<tensor, PatchField, Mesh>& gradD =
            mesh().lookupObject
            <
                GeometricField<tensor, PatchField, Mesh>
            >("grad(D)" + ext);

        epsilon = symm(gradD);
    }
    replacePlaneStress(epsilon, coeff, sigma);
}


template<template<class> class PatchField, class Mesh, class CoeffType>
void mechanicalLaw::updateEpsilon
(
    GeometricField<symmTensor, PatchField, Mesh>& epsilon,
    const tmp<CoeffType>& coeff,
    const GeometricField<symmTensor, PatchField, Mesh>& sigma
) const
{
    updateEpsilon(epsilon, coeff(), sigma);
}


template<template<class> class PatchField, class Mesh, class CoeffType>
void mechanicalLaw::updateEpsilon
(
    GeometricField<symmTensor, PatchField, Mesh>& epsilon,
    const tmp<CoeffType>& coeff,
    const GeometricField<symmTensor, PatchField, Mesh>& sigma,
    const GeometricField<symmTensor, PatchField, Mesh>& epsilonP
) const
{
    updateEpsilon(epsilon, coeff(), sigma, epsilonP);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
