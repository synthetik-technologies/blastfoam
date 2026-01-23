/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2024
     \\/     M anipulation  | Synthetik Applied Technology
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "StandardReconstructionScheme.H"
#include "fvc.H"


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
const Foam::surfaceScalarField&
Foam::StandardReconstructionScheme<Type>::lookupOrConstruct
(
    const word& fieldName,
    const scalar& value
) const
{
    if (!this->mesh_.template foundObject<surfaceScalarField>(fieldName))
    {
        surfaceScalarField* fPtr
        (
            new surfaceScalarField
            (
                IOobject
                (
                    fieldName,
                    this->mesh_.time().name(),
                    this->mesh_
                ),
                this->mesh_,
                value
            )
        );
        fPtr->store(fPtr);
    }

    // Reset value that may be flipped when balancing
    surfaceScalarField& f =
        this->mesh_.template lookupObjectRef<surfaceScalarField>(fieldName);
    f = value;

    return f;
}



template<class Type>
Foam::StandardReconstructionScheme<Type>::StandardReconstructionScheme
(
    const fvSchemes& schemes,
    const GeometricField<Type, fvPatchField, volMesh>& phi,
    Istream& is,
    const bool overwrite
)
:
    ReconstructionScheme<Type>(schemes, phi, is, overwrite),
    name_(is),
    ownNei_(ownNeiSurfaceFields::New(phi.mesh())),
    ownScheme_
    (
        surfaceInterpolationScheme<Type>::New
        (
            this->phi_.mesh(),
            ownNei_.own(),
            this->schemes_.interpolation(name_)
        )
    ),
    neiScheme_
    (
        surfaceInterpolationScheme<Type>::New
        (
            this->phi_.mesh(),
            ownNei_.nei(),
            this->schemes_.interpolation(name_)
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::StandardReconstructionScheme<Type>::~StandardReconstructionScheme()
{}

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::StandardReconstructionScheme<Type>::interpolateOwn() const
{
    return ownScheme_().interpolate(this->phi_);
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::StandardReconstructionScheme<Type>::interpolateNei() const
{
    return neiScheme_().interpolate(this->phi_);
}


// ************************************************************************* //
