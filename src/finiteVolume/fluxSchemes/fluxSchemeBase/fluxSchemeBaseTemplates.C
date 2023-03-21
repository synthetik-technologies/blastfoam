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
    const GeometricField<Type, fvPatchField, volMesh>& f
) const
{
    return interpolate<Type>(f, f.name());
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

    forAll(fOwn, facei)
    {
        Type& fi = ff[facei];
        const Type& fiOwn = fOwn[facei];
        const Type& fiNei = fNei[facei];
        for (label i = 0; i < nCmpts; i++)
        {
            setComponent(fi, i) = this->interpolate
            (
                component(fiOwn, i), component(fiNei, i),
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
                    this->interpolate
                    (
                        component(fiOwn, i), component(fiNei, i),
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
    const word& name
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;

    autoPtr<ReconstructionScheme<Type>> fLimiter
    (
        ReconstructionScheme<Type>::New(f, name)
    );

    tmp<fieldType> tfOwn, tfNei;
    fLimiter->interpolateOwnNei(tfOwn, tfNei);

    return interpolate(tfOwn(), tfNei(), name);
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> fluxSchemeBase::phaseInterpolate
(
    const volScalarField& alpha,
    const GeometricField<Type, fvPatchField, volMesh>& f,
    const scalar rAlpha
) const
{
    return phaseInterpolate<Type>(alpha, f, f.name(), rAlpha);
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> fluxSchemeBase::phaseInterpolate
(
    const volScalarField& alpha,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& fOwn,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& fNei,
    const word& name,
    const scalar rAlpha
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
    const labelList& owner = alpha.mesh().owner();
    const labelList& neighbour = alpha.mesh().neighbour();

    forAll(fOwn, facei)
    {
        Type& fi = ff[facei];
        const Type& fiOwn = fOwn[facei];
        const Type& fiNei = fNei[facei];
        const bool validOwn = alpha[owner[facei]] > rAlpha;
        const bool validNei = alpha[neighbour[facei]] > rAlpha;
        if (validOwn && validNei)
        {
            for (label i = 0; i < nCmpts; i++)
            {
                setComponent(fi, i) = this->interpolate
                (
                    component(fiOwn, i), component(fiNei, i),
                    facei
                );
            }
        }
        else if (validNei)
        {
            fi = fNei;
        }
        else
        {
            fi = fOwn;
        }

    }

    forAll(fOwn.boundaryField(), patchi)
    {
        Field<Type>& pff = ff.boundaryFieldRef()[patchi];
        const Field<Type>& pfOwn = fOwn.boundaryField()[patchi];
        const Field<Type>& pfNei = fNei.boundaryField()[patchi];
        if (alpha.boundaryField()[patchi].coupled())
        {
            const scalarField alphaOwn(alpha.boundaryField()[patchi].patchInternalField());
            const scalarField alphaNei(alpha.boundaryField()[patchi].patchNeighbourField());

            forAll(pff, facei)
            {
                Type& fi = pff[facei];
                const Type& fiOwn = pfOwn[facei];
                const Type& fiNei = pfNei[facei];
                const bool validOwn = alphaOwn[facei] > rAlpha;
                const bool validNei = alphaNei[facei] > rAlpha;

                if (validOwn && validNei)
                {
                    for (label i = 0; i < nCmpts; i++)
                    {
                        setComponent(fi, i) =
                            this->interpolate
                            (
                                component(fiOwn, i), component(fiNei, i),
                                facei, patchi
                            );
                    }
                }
                else if (validNei)
                {
                    fi = fiNei;
                }
                else
                {
                    fi = fiOwn;
                }
            }
        }

        else
        {
            pff = pfOwn;
        }
    }
    return tmpf;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> fluxSchemeBase::phaseInterpolate
(
    const volScalarField& alpha,
    const GeometricField<Type, fvPatchField, volMesh>& f,
    const word& name,
    const scalar rAlpha
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;

    autoPtr<ReconstructionScheme<Type>> fLimiter
    (
        ReconstructionScheme<Type>::New(f, name)
    );

    tmp<fieldType> tfOwn, tfNei;
    fLimiter->interpolateOwnNei(tfOwn, tfNei);

    return phaseInterpolate(alpha, tfOwn(), tfNei(), name, rAlpha);
}


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

    tmp<fieldType> tfOwn, tfNei;
    autoPtr<ReconstructionScheme<Type>> fLimiter
    (
        ReconstructionScheme<Type>::New(f, f.member(), f.group(), overwrite)
    );
    fLimiter->interpolateOwnNei(tfOwn, tfNei);

    return flux(tfOwn(), tfNei(), phi);
}


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

    tmp<surfaceScalarField> talphaRhoOwn, talphaRhoNei;
    autoPtr<ReconstructionScheme<scalar>> alphaRhoLimiter
    (
        ReconstructionScheme<scalar>::New
        (
            alphaRho,
            alphaRho.member(),
            alphaRho.group(),
            overwrite
        )
    );
    alphaRhoLimiter->interpolateOwnNei(talphaRhoOwn, talphaRhoNei);

    tmp<fieldType> tfOwn, tfNei;
    autoPtr<ReconstructionScheme<Type>> fLimiter
    (
        ReconstructionScheme<Type>::New(f, f.member(), f.group())
    );
    fLimiter->interpolateOwnNei(tfOwn, tfNei);
    tmp<fieldType> talphaRhofOwn(tfOwn*talphaRhoOwn);
    tmp<fieldType> talphaRhofNei(tfNei*talphaRhoNei);

    return flux(talphaRhofOwn(), talphaRhofNei(), phi);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, fvsPatchField, surfaceMesh>>
Foam::fluxSchemeBase::flux
(
    const GeometricField<Type, fvPatchField, volMesh>& f,
    const surfaceScalarField& alphaRhoOwn,
    const surfaceScalarField& alphaRhoNei,
    const surfaceScalarField& phi
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;

    tmp<fieldType> tfOwn, tfNei;
    autoPtr<ReconstructionScheme<Type>> fLimiter
    (
        ReconstructionScheme<Type>::New(f, f.member(), f.group())
    );
    fLimiter->interpolateOwnNei(tfOwn, tfNei);

    tmp<fieldType> talphaRhofOwn(tfOwn*alphaRhoOwn);
    tmp<fieldType> talphaRhofNei(tfNei*alphaRhoNei);

    return flux(talphaRhofOwn(), talphaRhofNei(), phi);
}


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


template<class Type>
Foam::tmp<Foam::GeometricField<Type, fvsPatchField, surfaceMesh>>
Foam::fluxSchemeBase::phaseFlux
(
    const volScalarField& alpha,
    const GeometricField<Type, fvPatchField, volMesh>& f,
    const surfaceScalarField& phi,
    const scalar rAlpha,
    const bool overwrite
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;

    tmp<fieldType> tfOwn, tfNei;
    autoPtr<ReconstructionScheme<Type>> fLimiter
    (
        ReconstructionScheme<Type>::New(f, f.member(), f.group(), overwrite)
    );
    fLimiter->interpolateOwnNei(tfOwn, tfNei);

    return phaseFlux(alpha, tfOwn(), tfNei(), phi, rAlpha);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, fvsPatchField, surfaceMesh>>
Foam::fluxSchemeBase::phaseFlux
(
    const volScalarField& alpha,
    const GeometricField<Type, fvPatchField, volMesh>& f,
    const volScalarField& alphaRho,
    const surfaceScalarField& phi,
    const scalar rAlpha,
    const bool overwrite
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;

    tmp<surfaceScalarField> talphaRhoOwn, talphaRhoNei;
    autoPtr<ReconstructionScheme<scalar>> alphaRhoLimiter
    (
        ReconstructionScheme<scalar>::New
        (
            alphaRho,
            alphaRho.member(),
            alphaRho.group(),
            overwrite
        )
    );
    alphaRhoLimiter->interpolateOwnNei(talphaRhoOwn, talphaRhoNei);

    tmp<fieldType> tfOwn, tfNei;
    autoPtr<ReconstructionScheme<Type>> fLimiter
    (
        ReconstructionScheme<Type>::New(f, f.member(), f.group())
    );
    fLimiter->interpolateOwnNei(tfOwn, tfNei);
    tmp<fieldType> talphaRhofOwn(tfOwn*talphaRhoOwn);
    tmp<fieldType> talphaRhofNei(tfNei*talphaRhoNei);

    return phaseFlux(alpha, talphaRhofOwn(), talphaRhofNei(), phi, rAlpha);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, fvsPatchField, surfaceMesh>>
Foam::fluxSchemeBase::phaseFlux
(
    const volScalarField& alpha,
    const GeometricField<Type, fvPatchField, volMesh>& f,
    const surfaceScalarField& alphaRhoOwn,
    const surfaceScalarField& alphaRhoNei,
    const surfaceScalarField& phi,
    const scalar rAlpha
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;

    tmp<fieldType> tfOwn, tfNei;
    autoPtr<ReconstructionScheme<Type>> fLimiter
    (
        ReconstructionScheme<Type>::New(f, f.member(), f.group())
    );
    fLimiter->interpolateOwnNei(tfOwn, tfNei);

    tmp<fieldType> talphaRhofOwn(tfOwn*alphaRhoOwn);
    tmp<fieldType> talphaRhofNei(tfNei*alphaRhoNei);

    return phaseFlux(alpha, talphaRhofOwn(), talphaRhofNei(), phi, rAlpha);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, fvsPatchField, surfaceMesh>>
Foam::fluxSchemeBase::phaseFlux
(
    const volScalarField& alpha,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& fOwn,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& fNei,
    const surfaceScalarField& phi,
    const scalar rAlpha
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
    const labelList& owner = alpha.mesh().owner();
    const labelList& neighbour = alpha.mesh().neighbour();

    forAll(fOwn, facei)
    {
        Type& fluxfi = fluxf[facei];
        const Type& fiOwn = fOwn[facei];
        const Type& fiNei = fNei[facei];
        const bool validOwn = alpha[owner[facei]] > rAlpha;
        const bool validNei = alpha[neighbour[facei]] > rAlpha;

        if (validOwn && validNei)
        {
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
        else if (validNei)
        {
            fluxfi = fiNei*phi[facei];
        }
        else
        {
            fluxfi = fiOwn*phi[facei];
        }
    }

    typename fieldType::Boundary& bfluxf = fluxf.boundaryFieldRef();
    forAll(fOwn.boundaryField(), patchi)
    {
        Field<Type>& pfluxf = bfluxf[patchi];
        const Field<Type>& pfOwn = fOwn.boundaryField()[patchi];
        const Field<Type>& pfNei = fNei.boundaryField()[patchi];
        const Field<scalar>& pphi = phi.boundaryField()[patchi];

        if (alpha.boundaryField()[patchi].coupled())
        {
            const scalarField alphaOwn(alpha.boundaryField()[patchi].patchInternalField());
            const scalarField alphaNei(alpha.boundaryField()[patchi].patchNeighbourField());

            forAll(pfluxf, facei)
            {
                Type& fluxfi = pfluxf[facei];
                const Type& fiOwn = pfOwn[facei];
                const Type& fiNei = pfNei[facei];
                const scalar& phii = pphi[facei];
                const bool validOwn = alphaOwn[facei] > rAlpha;
                const bool validNei = alphaNei[facei] > rAlpha;

                if (validOwn && validNei)
                {
                    for (label i = 0; i < nCmpts; i++)
                    {
                        setComponent(fluxfi, i) = this->calculateFlux
                        (
                            component(fiOwn, i), component(fiNei, i),
                            phii,
                            facei, patchi
                        );
                    }
                }
                else if (validNei)
                {
                    fluxfi = phii*fiNei;
                }
                else
                {
                    fluxfi = phii*fiOwn;
                }
            }
        }
        else
        {
            pfluxf = pphi*pfOwn;
        }
    }
    return tfluxf;
}


template<class Type>
void Foam::fluxSchemeBase::correctPhaseFields
(
    const volScalarField& alpha,
    GeometricField<Type, fvsPatchField, surfaceMesh>& fOwn,
    GeometricField<Type, fvsPatchField, surfaceMesh>& fNei,
    const scalar rAlpha
)
{
    const labelList& owner = alpha.mesh().owner();
    const labelList& neighbour = alpha.mesh().neighbour();

    forAll(fOwn, facei)
    {
        const bool validOwn = alpha[owner[facei]] > rAlpha;
        const bool validNei = alpha[neighbour[facei]] > rAlpha;

        if (validOwn && validNei)
        {
            // Do nothing
        }
        else if (validNei)
        {
            fOwn[facei] = fNei[facei];
        }
        else if (validOwn)
        {
            fNei[facei] = fOwn[facei];
        }
    }

    typename GeometricField<Type, fvsPatchField, surfaceMesh>::Boundary& bfOwn =
        fOwn.boundaryFieldRef();
    typename GeometricField<Type, fvsPatchField, surfaceMesh>::Boundary& bfNei =
        fNei.boundaryFieldRef();
    forAll(bfOwn, patchi)
    {
        Field<Type>& pfOwn = bfOwn[patchi];
        Field<Type>& pfNei = bfNei[patchi];
        const fvPatchField<scalar>& palpha = alpha.boundaryField()[patchi];

        if (alpha.boundaryField()[patchi].coupled())
        {
            const scalarField alphaOwn(palpha.patchInternalField());
            const scalarField alphaNei(palpha.patchNeighbourField());

            forAll(pfOwn, facei)
            {
                const bool validOwn = alphaOwn[facei] > rAlpha;
                const bool validNei = alphaNei[facei] > rAlpha;

                if (validOwn && validNei)
                {
                    // Do nothing
                }
                else if (validNei)
                {
                    pfOwn[facei] = pfNei[facei];
                }
                else if (validOwn)
                {
                    pfNei[facei] = pfOwn[facei];
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
