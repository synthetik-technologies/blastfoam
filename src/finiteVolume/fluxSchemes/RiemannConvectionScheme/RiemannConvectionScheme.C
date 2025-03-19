/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "RiemannConvectionScheme.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMatrices.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
RiemannConvectionScheme<Type>::interpolate
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return fluxSchemePtr_->interpolate(vf, vf.name());
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
RiemannConvectionScheme<Type>::flux
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return interpolate(faceFlux, vf)*faceFlux;
    if (&faceFlux == &(fluxSchemePtr_->phi()))
    {
        return fluxSchemePtr_->flux(vf, faceFlux, true);
    }
    else if
    (
        faceFlux.mesh().template foundObject<volScalarField>
        (
            IOobject::groupName("alphaRho", vf.group())
        )
    )
    {
        return fluxSchemePtr_->flux
        (
            vf,
            faceFlux.mesh().template lookupObject<volScalarField>
            (
                IOobject::groupName("alphaRho", vf.group())
            ),
            fluxSchemePtr_->phi(),
            false
        );
    }
    else if
    (
        faceFlux.mesh().template foundObject<volScalarField>
        (
            IOobject::groupName("rho", vf.group())
        )
    )
    {
        return fluxSchemePtr_->flux
        (
            vf,
            faceFlux.mesh().template lookupObject<volScalarField>
            (
                IOobject::groupName("rho", vf.group())
            ),
            fluxSchemePtr_->phi(),
            false
        );
    }
    else
    {
        return faceFlux*interpolate(faceFlux, vf);
    }
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
RiemannConvectionScheme<Type>::flux
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const volScalarField& alphaRho
) const
{
    return fluxSchemePtr_->flux(vf, alphaRho, faceFlux, true);
}


template<class Type>
tmp<fvMatrix<Type>>
RiemannConvectionScheme<Type>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{

    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            faceFlux.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    fvm += fvc::surfaceIntegrate(flux(faceFlux, vf));

    return tfvm;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
RiemannConvectionScheme<Type>::fvcDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> tConvection
    (
        fvc::surfaceIntegrate(flux(faceFlux, vf))
    );

    tConvection.ref().rename
    (
        "convection(" + faceFlux.name() + ',' + vf.name() + ')'
    );

    return tConvection;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
