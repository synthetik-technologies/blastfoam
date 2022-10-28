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

#include "fluxSchemeBase.H"
#include "ReconstructionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> fluxSchemeBase::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& f,
    const bool overwrite
) const
{
    return interpolate<Type>(f, f.name(), overwrite);
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> fluxSchemeBase::interpolate
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& fOwn,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& fNei,
    const word& name
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;
    label nCmpts = pTraits<Type>::nComponents;

    tmp<fieldType> tmpf
    (
        fieldType::New
        (
            name + "f",
            mesh_,
            dimensioned<Type>("0", fOwn.dimensions(), pTraits<Type>::zero)
        )
    );
    fieldType& ff = tmpf.ref();
    const bool isDensity = (fOwn.dimensions() == dimDensity);

    forAll(fOwn, facei)
    {
        Type& fi = ff[facei];
        const Type& fiOwn = fOwn[facei];
        const Type& fiNei = fNei[facei];
        for (label i = 0; i < nCmpts; i++)
        {
            setComponent(fi, i) = interpolate
            (
                component(fiOwn, i),
                component(fiNei, i),
                isDensity,
                facei
            );
        }
    }

    forAll(fOwn.boundaryField(), patchi)
    {
        Field<Type>& pff = ff.boundaryFieldRef()[patchi];
        const Field<Type>& pfOwn = fOwn.boundaryField()[patchi];
        const Field<Type>& pfNei = fNei.boundaryField()[patchi];

        forAll(pff, facei)
        {
            Type& fi = pff[facei];
            const Type& fiOwn = pfOwn[facei];
            const Type& fiNei = pfNei[facei];
            for (label i = 0; i < nCmpts; i++)
            {
                setComponent(fi, i) =
                    interpolate
                    (
                        component(fiOwn, i),
                        component(fiNei, i),
                        isDensity,
                        facei, patchi
                    );
            }
        }
    }
    return tmpf;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> fluxSchemeBase::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& f,
    const word& name,
    const bool overwrite
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;

    autoPtr<ReconstructionScheme<Type>> fLimiter
    (
        ReconstructionScheme<Type>::New(f, name)
    );

    tmp<fieldType> tfOwn;
    tmp<fieldType> tfNei;
    fLimiter->interpolateOwnNei(tfOwn, tfNei, overwrite);

    return interpolate(tfOwn(), tfNei(), name);
}


//- Own/nei interpolation of f
template<class Type>
Foam::tmp<Foam::GeometricField<Type, fvsPatchField, surfaceMesh>>
Foam::fluxSchemeBase::flux
(
    const GeometricField<Type, fvPatchField, volMesh>& f,
    const surfaceScalarField& phi,
    const bool overwrite
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;

    tmp<fieldType> tfOwn;
    tmp<fieldType> tfNei;
    autoPtr<ReconstructionScheme<Type>> fLimiter
    (
        ReconstructionScheme<Type>::New(f, f.member(), f.group())
    );
    fLimiter->interpolateOwnNei(tfOwn, tfNei, overwrite);

    return flux(tfOwn(), tfNei(), phi);
}

//- Own/nei interpolation of f and alphaRho
template<class Type>
Foam::tmp<Foam::GeometricField<Type, fvsPatchField, surfaceMesh>>
Foam::fluxSchemeBase::flux
(
    const GeometricField<Type, fvPatchField, volMesh>& f,
    const volScalarField& alphaRho,
    const surfaceScalarField& phi,
    const bool overwrite
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;

    tmp<surfaceScalarField> talphaRhoOwn;
    tmp<surfaceScalarField> talphaRhoNei;
    autoPtr<ReconstructionScheme<scalar>> alphaRhoLimiter
    (
        ReconstructionScheme<scalar>::New
        (
            alphaRho,
            alphaRho.member(),
            alphaRho.group()
        )
    );
    alphaRhoLimiter->interpolateOwnNei(talphaRhoOwn, talphaRhoNei, overwrite);

    tmp<fieldType> tfOwn;
    tmp<fieldType> tfNei;
    autoPtr<ReconstructionScheme<Type>> fLimiter
    (
        ReconstructionScheme<Type>::New(f, f.member(), f.group())
    );
    fLimiter->interpolateOwnNei(tfOwn, tfNei, overwrite);
    tmp<fieldType> talphaRhofOwn(tfOwn*talphaRhoOwn);
    tmp<fieldType> talphaRhofNei(tfNei*talphaRhoNei);

    return flux(talphaRhofOwn(), talphaRhofNei(), phi);
}

//- Own/nei interpolation of f
template<class Type>
Foam::tmp<Foam::GeometricField<Type, fvsPatchField, surfaceMesh>>
Foam::fluxSchemeBase::flux
(
    const GeometricField<Type, fvPatchField, volMesh>& f,
    const surfaceScalarField& alphaRhoOwn,
    const surfaceScalarField& alphaRhoNei,
    const surfaceScalarField& phi,
    const bool overwrite
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;

    tmp<fieldType> tfOwn;
    tmp<fieldType> tfNei;
    autoPtr<ReconstructionScheme<Type>> fLimiter
    (
        ReconstructionScheme<Type>::New(f, f.member(), f.group())
    );
    fLimiter->interpolateOwnNei(tfOwn, tfNei);

    tmp<fieldType> talphaRhofOwn(tfOwn*alphaRhoOwn);
    tmp<fieldType> talphaRhofNei(tfNei*alphaRhoNei);

    return flux(talphaRhofOwn(), talphaRhofNei(), phi);
}


//- Own/nei interpolation of f
template<class Type>
Foam::tmp<Foam::GeometricField<Type, fvsPatchField, surfaceMesh>>
Foam::fluxSchemeBase::flux
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& fOwn,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& fNei,
    const surfaceScalarField& phi
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;
    label nCmpts = pTraits<Type>::nComponents;

    tmp<fieldType> tfluxf
    (
        fieldType::New
        (
            "flux(" + fOwn.name() + "|" + fNei.name() + "," + phi.name() + ")",
            mesh_,
            dimensioned<Type>("0", fOwn.dimensions()*phi.dimensions(), Zero)
        )
    );
    fieldType& fluxf = tfluxf.ref();

    forAll(fOwn, facei)
    {
        Type& fluxfi = fluxf[facei];
        const Type& fiOwn = fOwn[facei];
        const Type& fiNei = fNei[facei];
        for (label i = 0; i < nCmpts; i++)
        {
            setComponent(fluxfi, i) = calculateFlux
            (
                component(fiOwn, i), component(fiNei, i),
                phi[facei],
                facei
            );
        }
    }

    typename fieldType::Boundary& bfluxf = fluxf.boundaryFieldRef();
    forAll(fOwn.boundaryField(), patchi)
    {
        Field<Type>& pfluxf = bfluxf[patchi];
        const Field<Type>& pfOwn = fOwn.boundaryField()[patchi];
        const Field<Type>& pfNei = fNei.boundaryField()[patchi];
        const Field<scalar>& pphi = phi.boundaryField()[patchi];

        forAll(pfluxf, facei)
        {
            Type& fluxfi = pfluxf[facei];
            const Type& fiOwn = pfOwn[facei];
            const Type& fiNei = pfNei[facei];
            const scalar& phii = pphi[facei];
            for (label i = 0; i < nCmpts; i++)
            {
                setComponent(fluxfi, i) = calculateFlux
                (
                    component(fiOwn, i), component(fiNei, i),
                    phii,
                    facei, patchi
                );
            }
        }
    }
    return tfluxf;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
