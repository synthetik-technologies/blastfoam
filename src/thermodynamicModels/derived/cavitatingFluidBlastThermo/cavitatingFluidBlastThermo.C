/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a cavitating material
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

#include "cavitatingFluidBlastThermo.H"
#include "fvc.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
void Foam::cavitatingFluidBlastThermo<Thermo>::calculate()
{
    const typename Thermo::thermoType1& tl(*this);
    const typename Thermo::thermoType2& tv(*this);

    volScalarField::Internal& CpI = this->CpRef().internalFieldRef();
    volScalarField::Internal& CvI = this->CvRef().internalFieldRef();
    volScalarField::Internal& muI = this->muRef().internalFieldRef();
    volScalarField::Internal& kappaI = this->kappaRef().internalFieldRef();
    volScalarField::Internal& pI = this->pRef().internalFieldRef();
    volScalarField::Internal& cI = this->speedOfSoundRef().internalFieldRef();

    forAll(this->rho_, celli)
    {
        scalar& ei(this->heRef()[celli]);
        scalar& Ti(this->TRef()[celli]);


        const scalar rhoi(this->rho_[celli]);
        const scalar xv = this->cellx(celli);
        const scalar xl = 1.0 - xv;

        if (xl < small)
        {
            Ti = tv.TRhoE(Ti, rhoi, ei);
            if (Ti < this->TLow_)
            {
                ei = tv.Es(rhoi, ei, this->TLow_);
                Ti = this->TLow_;
            }

            const scalar pi = tv.p(rhoi, ei, Ti);
            pI[celli] = pi;
            CpI[celli] = tv.Cp(rhoi, ei, Ti);
            CvI[celli] = tv.Cv(rhoi, ei, Ti);
            muI[celli] = tv.mu(rhoi, ei, Ti);
            kappaI[celli] = tv.kappa(rhoi, ei, Ti);
            cI[celli] = sqrt(max(tv.cSqr(pi, rhoi, ei, Ti), small));
        }
        else if (xv < small)
        {
            Ti = tl.TRhoE(Ti, rhoi, ei);
            if (Ti < this->TLow_)
            {
                ei = tl.Es(rhoi, ei, this->TLow_);
                Ti = this->TLow_;
            }

            const scalar pi = tl.p(rhoi, ei, Ti);

            pI[celli] = pi;
            CpI[celli] = tl.Cp(rhoi, ei, Ti);
            CvI[celli] = tl.Cv(rhoi, ei, Ti);
            muI[celli] = tl.mu(rhoi, ei, Ti);
            kappaI[celli] = tl.kappa(rhoi, ei, Ti);
            cI[celli] = sqrt(max(tl.cSqr(pi, rhoi, ei, Ti), small));
        }
        else
        {
            const scalar pSat = pSat_->pSat(Ti);
            const scalar rhoSatv = rhoSatv_.lookup(Ti);
            const scalar rhoSatl = rhoSatl_.lookup(Ti);

            const scalar fv = rhoSatv/rhoi;
            const scalar fl = rhoSatl/rhoi;

            Ti =
                tv.TRhoE(Ti, rhoi, ei)*xv
              + tl.TRhoE(Ti, rhoi, ei)*xl;
            if (Ti < this->TLow_)
            {
                ei =
                    tv.Es(rhoi, ei, this->TLow_)*xv*fv
                  + tl.Es(rhoi, ei, this->TLow_)*xl*fl;
                Ti = this->TLow_;
            }

            const scalar pi = pSat;

            pI[celli] = pi;
            CpI[celli] =
                tv.Cp(rhoi, ei, Ti)*xv*fv
              + tl.Cp(rhoi, ei, Ti)*xl*fl;
            CvI[celli] =
                tv.Cv(rhoi, ei, Ti)*xv*fv
              + tl.Cv(rhoi, ei, Ti)*xl*fl;
            muI[celli] =
                tv.mu(rhoi, ei, Ti)*xv
              + tl.mu(rhoi, ei, Ti)*xl;
            kappaI[celli] =
                tv.kappa(rhoi, ei, Ti)*xv
              + tl.kappa(rhoi, ei, Ti)*xl;
            cI[celli] =
                sqrt
                (
                    1.0
                   /(
                        xv/(rhoSatv*max(tv.cSqr(pi, rhoi, ei, Ti), small))
                      + xl/(rhoSatl*max(tl.cSqr(pi, rhoi, ei, Ti), small))
                    )/rhoi
                );
        }
    }

    this->pRef().correctBoundaryConditions();

    volScalarField::Boundary& bhe = this->heRef().boundaryFieldRef();
    volScalarField::Boundary& bT = this->TRef().boundaryFieldRef();
    volScalarField::Boundary& bCp = this->CpRef().boundaryFieldRef();
    volScalarField::Boundary& bCv = this->CvRef().boundaryFieldRef();
    volScalarField::Boundary& bmu = this->muRef().boundaryFieldRef();
    volScalarField::Boundary& bkappa = this->kappaRef().boundaryFieldRef();
    volScalarField::Boundary& bc = this->speedOfSoundRef().boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pp =
            this->pRef().boundaryField()[patchi];
        const fvPatchScalarField& px = x_.boundaryField()[patchi];

        fvPatchScalarField& pT = bT[patchi];
        fvPatchScalarField& phe = bhe[patchi];
        fvPatchScalarField& pCp = bCp[patchi];
        fvPatchScalarField& pCv = bCv[patchi];
        fvPatchScalarField& pmu = bmu[patchi];
        fvPatchScalarField& pkappa = bkappa[patchi];
        fvPatchScalarField& pc = bc[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const scalar rhoi(prho[facei]);
                scalar& ei = phe[facei];
                scalar& Ti = pT[facei];

                const scalar xv = px[facei];
                const scalar xl = 1.0 - xv;

                const scalar pi(pp[facei]);

                if (xl < this->residualFac_)
                {
                    ei = tv.Es(rhoi, ei, Ti);
                    pCp[facei] = tv.Cp(rhoi, ei, Ti);
                    pCv[facei] = tv.Cv(rhoi, ei, Ti);
                    pmu[facei] = tv.mu(rhoi, ei, Ti);
                    pkappa[facei] = tv.kappa(rhoi, ei, Ti);
                    pc[facei] =
                        sqrt(max(tv.cSqr(pi, rhoi, ei, Ti), small));
                }
                else if (xv < this->residualFac_)
                {
                    ei = tl.Es(rhoi, ei, Ti);
                    pCp[facei] = tl.Cp(rhoi, ei, Ti);
                    pCv[facei] = tl.Cv(rhoi, ei, Ti);
                    pmu[facei] = tl.mu(rhoi, ei, Ti);
                    pkappa[facei] = tl.kappa(rhoi, ei, Ti);
                    pc[facei] =
                        sqrt(max(tl.cSqr(pi, rhoi, ei, Ti), small));
                }
                else
                {
                    const scalar rhoSatv = rhoSatv_.lookup(Ti);
                    const scalar rhoSatl = rhoSatl_.lookup(Ti);

                    const scalar fv = rhoSatv/rhoi;
                    const scalar fl = rhoSatl/rhoi;

                    ei =
                        tv.Es(rhoi, ei, Ti)*xv
                      + tl.Es(rhoi, ei , Ti)*xl;
                    pCp[facei] =
                        tv.Cp(rhoi, ei, Ti)*xv*fv
                      + tl.Cp(rhoi, ei, Ti)*xl*fl;
                    pCv[facei] =
                        tv.Cv(rhoi, ei, Ti)*xv*fv
                      + tl.Cv(rhoi, ei, Ti)*xl*fl;
                    pmu[facei] =
                        tv.mu(rhoi, ei, Ti)*xv
                      + tl.mu(rhoi, ei, Ti)*xl;
                    pkappa[facei] =
                        tv.kappa(rhoi, ei, Ti)*xv
                      + tl.kappa(rhoi, ei, Ti)*xl;
                    pc[facei] =
                        sqrt
                        (
                            1.0
                           /(
                                xv
                                /(
                                    rhoSatv
                                    *max(tv.cSqr(pi, rhoi, ei, Ti), small)
                                )
                              + xl
                               /(
                                    rhoSatl
                                   *max(tl.cSqr(pi, rhoi, ei, Ti), small)
                                )
                            )/rhoi
                        );
                }
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const scalar rhoi(prho[facei]);

                const scalar xv = px[facei];
                const scalar xl = 1.0 - xv;

                const scalar pi(pp[facei]);

                scalar& ei = phe[facei];
                scalar& Ti = pT[facei];

                if (xl < this->residualFac_)
                {
                    Ti = tv.TRhoE(Ti, rhoi, ei);
                    if (Ti < this->TLow_)
                    {
                        Ti = this->TLow_;
                        ei = tv.Es(rhoi, ei, Ti);
                    }
                    pCp[facei] = tv.Cp(rhoi, ei, Ti);
                    pCv[facei] = tv.Cv(rhoi, ei, Ti);
                    pmu[facei] = tv.mu(rhoi, ei, Ti);
                    pkappa[facei] = tv.kappa(rhoi, ei, Ti);
                    pc[facei] =
                        sqrt(max(tv.cSqr(pi, rhoi, ei, Ti), small));
                }
                else if (xv < this->residualFac_)
                {
                    Ti = tl.TRhoE(Ti, rhoi, ei);
                    if (Ti < this->TLow_)
                    {
                        Ti = this->TLow_;
                        ei = tl.Es(rhoi, ei, Ti);
                    }
                    pCp[facei] = tl.Cp(rhoi, ei, Ti);
                    pCv[facei] = tl.Cv(rhoi, ei, Ti);
                    pmu[facei] = tl.mu(rhoi, ei, Ti);
                    pkappa[facei] = tl.kappa(rhoi, ei, Ti);
                    pc[facei] =
                        sqrt(max(tl.cSqr(pi, rhoi, ei, Ti), small));
                }
                else
                {
                    const scalar rhoSatv = rhoSatv_.lookup(Ti);
                    const scalar rhoSatl = rhoSatl_.lookup(Ti);

                    const scalar fv = rhoSatv/rhoi;
                    const scalar fl = rhoSatl/rhoi;

                    Ti =
                        tv.TRhoE(Ti, rhoi, ei)*xv
                      + tl.TRhoE(Ti, rhoi, ei)*xl;
                    if (Ti < this->TLow_)
                    {
                        Ti = this->TLow_;
                        ei =
                            tv.Es(rhoi, ei, Ti)*xv
                          + tl.Es(rhoi, ei, Ti)*xl;
                    }
                    pCp[facei] =
                        tv.Cp(rhoi, ei, Ti)*xv*fv
                      + tl.Cp(rhoi, ei, Ti)*xl*fl;
                    pCv[facei] =
                        tv.Cv(rhoi, ei, Ti)*xv*fv
                      + tl.Cv(rhoi, ei, Ti)*xl*fl;
                    pmu[facei] =
                        tv.mu(rhoi, ei, Ti)*xv
                      + tl.mu(rhoi, ei, Ti)*xl;
                    pkappa[facei] =
                        tv.kappa(rhoi, ei, Ti)*xv
                      + tl.kappa(rhoi, ei, Ti)*xl;
                    pc[facei] =
                        sqrt
                        (
                            1.0
                           /(
                                xv
                               /(
                                    rhoSatv
                                   *max(tv.cSqr(pi, rhoi, ei, Ti), small)
                                )
                              + xl
                               /(
                                    rhoSatl
                                   *max(tl.cSqr(pi, rhoi, ei, Ti), small)
                                )
                            )/rhoi
                        );
                }
            }
        }

    }
}



template<class Thermo>
void Foam::cavitatingFluidBlastThermo<Thermo>::calculate
(
    const volScalarField& alpha,
    const volScalarField& he,
    const volScalarField& T,
    volScalarField& alphaCp,
    volScalarField& alphaCv,
    volScalarField& alphaMu,
    volScalarField& alphaAlphah,
    volScalarField& pXiSum,
    volScalarField& XiSum
)
{
    const typename Thermo::thermoType1& tl(*this);
    const typename Thermo::thermoType2& tv(*this);

    forAll(alpha, celli)
    {
        const scalar xv = this->cellx(celli);
        const scalar xl = 1.0 - xv;
        const scalar alphai = alpha[celli];
        const scalar rhoi = this->rho_[celli];
        const scalar ei = he[celli];
        const scalar Ti = T[celli];
        if
        (
            alphai > this->residualAlpha_.value()
         && rhoi > this->residualRho_.value()
        )
        {
            scalar Gamma;
            scalar pi;

            if (xl < this->residualFac_)
            {
                alphaCp[celli] += tv.Cp(rhoi, ei, Ti)*alphai;
                alphaCv[celli] += tv.Cv(rhoi, ei, Ti)*alphai;
                alphaMu[celli] += tv.mu(rhoi, ei, Ti)*alphai;
                alphaAlphah[celli] +=
                    tv.kappa(rhoi, ei, Ti)/tv.Cp(rhoi, ei, Ti)*alphai;
                Gamma = tv.Gamma(rhoi, ei, Ti);
                pi = tv.p(rhoi, ei, Ti);
            }
            else if (xv < this->residualFac_)
            {
                alphaCp[celli] += tl.Cp(rhoi, ei, Ti)*alphai;
                alphaCv[celli] += tl.Cv(rhoi, ei, Ti)*alphai;
                alphaMu[celli] += tl.mu(rhoi, ei, Ti)*alphai;
                alphaAlphah[celli] +=
                    tl.kappa(rhoi, ei, Ti)/tl.Cp(rhoi, ei, Ti)*alphai;

                Gamma = tl.Gamma(rhoi, ei, Ti);
                pi = tl.p(rhoi, ei, Ti);
            }
            else
            {
                const scalar pSat = pSat_->pSat(Ti);
                const scalar rhoSatv = rhoSatv_.lookup(Ti);
                const scalar rhoSatl = rhoSatl_.lookup(Ti);

                const scalar fv = rhoSatv/max(rhoi, 1e-10);
                const scalar fl = rhoSatl/max(rhoi, 1e-10);

                alphaCp[celli] +=
                    (
                        tv.Cp(rhoi, ei, Ti)*xv*fv
                      + tl.Cp(rhoi, ei, Ti)*xl*fl
                    )*alphai;
                alphaCv[celli] +=
                    (
                        tv.Cv(rhoi, ei, Ti)*xv*fv
                      + tl.Cv(rhoi, ei, Ti)*xl*fl
                    )*alphai;
                alphaMu[celli] +=
                    (
                        tv.mu(rhoi, ei, Ti)*xv
                      + tl.mu(rhoi, ei, Ti)*xl
                    )*alphai;
                alphaAlphah[celli] +=
                    (
                        tv.kappa(rhoi, ei, Ti)/tv.Cp(rhoi, ei, Ti)*xv
                      + tl.kappa(rhoi, ei, Ti)/tl.Cp(rhoi, ei, Ti)*xl
                    )*alphai;

                Gamma =
                    tv.Gamma(rhoi, ei, Ti)*xv + tv.Gamma(rhoi, ei, Ti)*xl;
                pi = pSat;
            }
            scalar Xii = alphai/(Gamma - 1.0);
            pXiSum[celli] += pi*Xii;
            XiSum[celli] += Xii;
        }
    }

    forAll(alpha.boundaryField(), patchi)
    {
        const fvPatchScalarField& palpha = alpha.boundaryField()[patchi];
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pT = T.boundaryField()[patchi];
        const fvPatchScalarField& phe = he.boundaryField()[patchi];
        const fvPatchScalarField& px = x_.boundaryField()[patchi];

        fvPatchScalarField& palphaCp = alphaCp.boundaryFieldRef()[patchi];
        fvPatchScalarField& palphaCv = alphaCv.boundaryFieldRef()[patchi];
        fvPatchScalarField& palphaMu = alphaMu.boundaryFieldRef()[patchi];
        fvPatchScalarField& palphaAlphah =
            alphaAlphah.boundaryFieldRef()[patchi];
        fvPatchScalarField& ppXiSum = pXiSum.boundaryFieldRef()[patchi];
        fvPatchScalarField& pxiSum = XiSum.boundaryFieldRef()[patchi];

        forAll(palpha, facei)
        {
            const scalar alphai = palpha[facei];
            const scalar rhoi = prho[facei];
            if
            (
                alphai > this->residualAlpha_.value()
             && rhoi > this->residualRho_.value()
            )
            {
                const scalar xv = px[facei];
                const scalar xl = 1.0 - xv;
                const scalar ei = phe[facei];
                const scalar Ti = pT[facei];
                scalar Gamma = 2.0;
                scalar pi = 0.0;

                if (xl < this->residualFac_)
                {
                    palphaCp[facei] += tv.Cp(rhoi, ei, Ti)*alphai;
                    palphaCv[facei] += tv.Cv(rhoi, ei, Ti)*alphai;
                    palphaMu[facei] += tv.mu(rhoi, ei, Ti)*alphai;
                    palphaAlphah[facei] +=
                        tv.kappa(rhoi, ei, Ti)/tv.Cp(rhoi, ei, Ti)*alphai;

                    Gamma = tv.Gamma(rhoi, ei, Ti);
                    pi = tv.p(rhoi, ei, Ti);
                }
                else if (xv < this->residualFac_)
                {
                    palphaCp[facei] += tl.Cp(rhoi, ei, Ti)*alphai;
                    palphaCv[facei] += tl.Cv(rhoi, ei, Ti)*alphai;
                    palphaMu[facei] += tl.mu(rhoi, ei, Ti)*alphai;
                    palphaAlphah[facei] +=
                        tl.kappa(rhoi, ei, Ti)/tl.Cp(rhoi, ei, Ti)*alphai;

                    Gamma = tl.Gamma(rhoi, ei, Ti);
                    pi = tl.p(rhoi, ei, Ti);
                }
                else
                {
                    const scalar pSat = pSat_->pSat(Ti);
                    const scalar rhoSatv = rhoSatv_.lookup(Ti);
                    const scalar rhoSatl = rhoSatl_.lookup(Ti);

                    const scalar fv = rhoSatv/max(rhoi, 1e-10);
                    const scalar fl = rhoSatl/max(rhoi, 1e-10);

                    palphaCp[facei] +=
                        (
                            tv.Cp(rhoi, ei, Ti)*xv*fv
                          + tl.Cp(rhoi, ei, Ti)*xl*fl
                        )*alphai;
                    palphaCv[facei] +=
                        (
                            tv.Cv(rhoi, ei, Ti)*xv*fv
                          + tl.Cv(rhoi, ei, Ti)*xl*fl
                        )*alphai;
                    palphaMu[facei] +=
                        (
                            tv.mu(rhoi, ei, Ti)*xv
                          + tl.mu(rhoi, ei, Ti)*xl
                        )*alphai;
                    palphaAlphah[facei] +=
                        (
                            tv.kappa(rhoi, ei, Ti)/tv.Cp(rhoi, ei, Ti)*xv
                          + tl.kappa(rhoi, ei, Ti)/tl.Cp(rhoi, ei, Ti)*xl
                        )*alphai;

                    Gamma =
                        tv.Gamma(rhoi, ei, Ti)*xv + tl.Gamma(rhoi, ei, Ti)*xl;
                    pi = pSat;
                }
                scalar Xii = alphai/(Gamma - 1.0);
                ppXiSum[facei] += pi*Xii;
                pxiSum[facei] += Xii;
            }
        }
    }
}


template<class Thermo>
void Foam::cavitatingFluidBlastThermo<Thermo>::calculateSpeedOfSound
(
    const volScalarField& alpha,
    volScalarField& cSqrRhoXiSum
)
{
    const typename Thermo::thermoType1& tl(*this);
    const typename Thermo::thermoType2& tv(*this);

    forAll(this->rho_, celli)
    {
        const scalar alphai = alpha[celli];
        const scalar rhoi = this->rho_[celli];
        if
        (
            alphai > this->residualAlpha_.value()
         && rhoi > this->residualRho_.value()
        )
        {
            const scalar xv = this->cellx(celli);
            const scalar xl = 1.0 - xv;
            const scalar pi = this->p_[celli];
            const scalar ei = this->e_[celli];
            const scalar Ti = this->T_[celli];
            scalar cSqr;
            scalar Gamma;

            if (xl < this->residualFac_)
            {
                cSqr = tv.cSqr(pi, rhoi, ei, Ti);
                Gamma = tv.Gamma(rhoi, ei, Ti);
            }
            else if (xv < this->residualFac_)
            {
                cSqr = tl.cSqr(pi, rhoi, ei, Ti);
                Gamma = tl.Gamma(rhoi, ei, Ti);
            }
            else
            {
                const scalar rhoSatv = rhoSatv_.lookup(Ti);
                const scalar rhoSatl = rhoSatl_.lookup(Ti);

                cSqr =
                    1.0
                   /(
                        xv/(rhoSatv*max(tv.cSqr(pi, rhoi, ei, Ti), small))
                      + xl/(rhoSatl*max(tl.cSqr(pi, rhoi, ei, Ti), small))
                    )/max(rhoi, 1e-10);
                Gamma =
                    tv.Gamma(rhoi, ei, Ti)*xv + tl.Gamma(rhoi, ei, Ti)*xl;
            }
            cSqrRhoXiSum[celli] += cSqr*rhoi*alphai/(Gamma - 1.0);
        }
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& palpha = alpha.boundaryField()[patchi];
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        const fvPatchScalarField& phe = this->e_.boundaryField()[patchi];
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& px = x_.boundaryField()[patchi];
        fvPatchScalarField& pcSqrRhoXiSum =
            cSqrRhoXiSum.boundaryFieldRef()[patchi];

        forAll(pT, facei)
        {
            const scalar alphai = palpha[facei];
            const scalar rhoi = prho[facei];
            if
            (
                alphai > this->residualAlpha_.value()
              && rhoi > this->residualRho_.value()
            )
            {
                const scalar xv = px[facei];
                const scalar xl = 1.0 - xv;
                const scalar pi = pp[facei];
                const scalar ei = phe[facei];
                const scalar Ti = pT[facei];
                scalar cSqr;
                scalar Gamma;

                if (xl < this->residualFac_)
                {
                    cSqr = tv.cSqr(pi, rhoi, ei, Ti);
                    Gamma = tv.Gamma(rhoi, ei, Ti);
                }
                else if (xv < this->residualFac_)
                {
                    cSqr = tl.cSqr(pi, rhoi, ei, Ti);
                    Gamma = tl.Gamma(rhoi, ei, Ti);
                }
                else
                {
                    const scalar rhoSatv = rhoSatv_.lookup(Ti);
                    const scalar rhoSatl = rhoSatl_.lookup(Ti);

                    cSqr =
                        1.0
                       /(
                            xv/(rhoSatv*max(tv.cSqr(pi, rhoi, ei, Ti), small))
                          + xl/(rhoSatl*max(tl.cSqr(pi, rhoi, ei, Ti), small))
                        )/max(rhoi, 1e-10);
                    Gamma =
                        tv.Gamma(rhoi, ei, Ti)*xv + tl.Gamma(rhoi, ei, Ti)*xl;
                }
                pcSqrRhoXiSum[facei] += cSqr*rhoi*alphai/(Gamma - 1.0);
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::cavitatingFluidBlastThermo<Thermo>::cavitatingFluidBlastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
:
    Thermo
    (
        mesh,
        dict,
        dict.subDict("liquid"),
        dict.subDict("vapor"),
        phaseName,
        masterName
    ),
    x_
    (
        IOobject
        (
            IOobject::groupName("x", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        0.0,
        "zeroGradient"
    ),
    pSat_(saturationPressureModel::New("pSat", dict)),
    hv_("Hv", dimEnergy/dimMass, dict.lookup<scalar>("Hv"))
{
    const typename Thermo::thermoType1& tl(*this);
    const typename Thermo::thermoType2& tv(*this);
    scalarList Ts(1001);
    scalarList rhol(1001);
    scalarList rhov(1001);
    forAll(Ts, i)
    {
        Ts[i] = 273.15 + scalar(i);
        const scalar pSat = pSat_->pSat(Ts[i]);
        rhov[i] = tv.rhoPTOffset(1.0, pSat, Ts[i], 0.0);
        rhol[i] = tl.rhoPT(1000.0, pSat, Ts[i]);
    }
    rhoSatl_.set
    (
        Ts,
        rhol,
        "none",
        "none",
        "linearClamp",
        false
    );
    rhoSatv_.set
    (
        Ts,
        rhov,
        "none",
        "none",
        "linearClamp",
        false
    );

    //- Initialize the density using the pressure and temperature
    //  This is only done at the first time step (Not on restart)
    if
    (
        max(this->rho_).value() <= 0
     || (
            Thermo::thermoType1::solid()
         && dict.lookupOrDefault<Switch>("calculateDensity", false)
         && !mesh.time().restart()
        )
    )
    {
        updateRho(Thermo::baseThermo::p());
    }
    this->initializeFields();
}


template<class Thermo>
void Foam::cavitatingFluidBlastThermo<Thermo>::initializeModels()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::cavitatingFluidBlastThermo<Thermo>::~cavitatingFluidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::cavitatingFluidBlastThermo<Thermo>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    calculate();

    if (debug)
    {
        Info<< "    Finished" << endl;
    }
}


template<class Thermo>
void Foam::cavitatingFluidBlastThermo<Thermo>::update()
{
    dxdt_.clear();
    alphaRhoXOld_.clear();
    deltaAlphaRhoX_.clear();

    const fvMesh& mesh = this->rho_.mesh();

    const volScalarField& alphaRho = mesh.lookupObject<volScalarField>
    (
        this->phaseName().empty()
      ? "rho"
      : IOobject::groupName("alphaRho", this->phaseName())
    );
    const surfaceScalarField& alphaRhoPhi = mesh.lookupObject<surfaceScalarField>
    (
        this->phaseName().empty()
      ? "rhoPhi"
      : IOobject::groupName("alphaRhoPhi", this->phaseName())
    );

    alphaRhoXOld_ = x_*alphaRho;
    volScalarField xOld(x_);
    this->storeAndBlendOld(xOld, false);
    volScalarField xNew(x_);
    forAll(x_, celli)
    {
        const scalar Ti = this->T_[celli];
        const scalar rhoi = this->rho_[celli];
        const scalar rhoSatv = rhoSatv_.lookup(Ti);
        const scalar rhoSatl = rhoSatl_.lookup(Ti);
        xNew[celli] = (rhoi - rhoSatl)/(rhoSatv - rhoSatl);
    }

    volScalarField deltax((xNew - x_)/mesh.time().deltaT());
    deltaAlphaRhoX_ =
        fvc::div(alphaRhoPhi, x_)
      - deltax*alphaRho;

    xNew.maxMin(0.0, 1.0);
    deltax = (xNew - xOld)/mesh.time().deltaT();
    dxdt_ = this->calcAndStoreDelta(deltax);

}


template<class Thermo>
void Foam::cavitatingFluidBlastThermo<Thermo>::solve()
{
    if (!alphaRhoXOld_.valid())
    {
        return;
    }

    const fvMesh& mesh = this->rho_.mesh();
    const volScalarField& alphaRho = mesh.lookupObject<volScalarField>
    (
        this->phaseName().empty()
      ? "rho"
      : IOobject::groupName("alphaRho", this->phaseName())
    );

    dimensionedScalar dT(mesh.time().deltaT());

    this->storeAndBlendOld(alphaRhoXOld_.ref());
    this->storeAndBlendDelta(deltaAlphaRhoX_.ref());

    //- Update lambda to include advection and reaction
    //  d(alpha rho lambda)/dt = alpha rho d(lambda)/dt + lambda d(alpha rho)/dt
    x_ =
        (alphaRhoXOld_ - deltaAlphaRhoX_*dT)
       /max(alphaRho, this->residualAlpha()*this->residualRho());
    x_.maxMin(0.0, 1.0);
    x_.correctBoundaryConditions();
}


template<class Thermo>
void Foam::cavitatingFluidBlastThermo<Thermo>::postUpdate()
{}


template<class Thermo>
void Foam::cavitatingFluidBlastThermo<Thermo>::updateRho(const volScalarField& p)
{
    const typename Thermo::thermoType1& tl(*this);
    const typename Thermo::thermoType2& tv(*this);

    scalarField& rhoI = this->rho_.primitiveFieldRef();
    forAll(this->rho_, celli)
    {
        const scalar xv = this->cellx(celli);
        const scalar xl = 1.0 - xv;

        if (xl < this->residualFac_)
        {
            rhoI[celli] = tv.rhoPT(rhoI[celli], p[celli], this->T_[celli]);
        }
        else if (xv < this->residualFac_)
        {
            rhoI[celli] = tl.rhoPT(rhoI[celli], p[celli], this->T_[celli]);
        }
        else
        {
            rhoI[celli] =
                tv.rhoPT(rhoI[celli], p[celli], this->T_[celli])*xv
              + tl.rhoPT(rhoI[celli], p[celli], this->T_[celli])*xl;
        }
    }

    volScalarField::Boundary& brho = this->rho_.boundaryFieldRef();

    forAll(brho, patchi)
    {
        scalarField& prho = brho[patchi];
        const scalarField& pT = this->T_.boundaryField()[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const fvPatchScalarField& px = x_.boundaryField()[patchi];

        forAll(prho, facei)
        {
            const scalar xv = px[facei];
            const scalar xl = 1.0 - xv;
            if (xl < this->residualFac_)
            {
                prho[facei] = tv.rhoPT(prho[facei], pp[facei], pT[facei]);
            }
            else if (xv < this->residualFac_)
            {
                prho[facei] = tl.rhoPT(prho[facei], pp[facei], pT[facei]);
            }
            else
            {
                prho[facei] =
                    tv.rhoPT(prho[facei], pp[facei], pT[facei])*xv
                  + tl.rhoPT(prho[facei], pp[facei], pT[facei])*xl;
            }
        }
    }
}


template<class Thermo>
void Foam::cavitatingFluidBlastThermo<Thermo>::updateRho
(
    const volScalarField& alpha,
    const volScalarField& p
)
{
    const typename Thermo::thermoType1& tl(*this);
    const typename Thermo::thermoType2& tv(*this);

    scalarField& rhoI = this->rho_.primitiveFieldRef();
    forAll(this->rho_, celli)
    {
        if (alpha[celli] > this->residualAlpha_.value())
        {
            const scalar xv = this->cellx(celli);
            const scalar xl = 1.0 - xv;

            if (xl < this->residualFac_)
            {
                rhoI[celli] = tv.rhoPT(rhoI[celli], p[celli], this->T_[celli]);
            }
            else if (xv < this->residualFac_)
            {
                rhoI[celli] = tl.rhoPT(rhoI[celli], p[celli], this->T_[celli]);
            }
            else
            {
                rhoI[celli] =
                    tv.rhoPT(rhoI[celli], p[celli], this->T_[celli])*xv
                  + tl.rhoPT(rhoI[celli], p[celli], this->T_[celli])*xl;
            }
        }
    }

    volScalarField::Boundary& brho = this->rho_.boundaryFieldRef();

    forAll(brho, patchi)
    {
        scalarField& prho = brho[patchi];
        const scalarField& palpha = alpha.boundaryField()[patchi];
        const scalarField& pT = this->T_.boundaryField()[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const fvPatchScalarField& px = x_.boundaryField()[patchi];

        forAll(prho, facei)
        {
            if (palpha[facei] > this->residualAlpha_.value())
            {
                const scalar xv = px[facei];
                const scalar xl = 1.0 - xv;
                if (xl < this->residualFac_)
                {
                    prho[facei] = tv.rhoPT(prho[facei], pp[facei], pT[facei]);
                }
                else if (xv < this->residualFac_)
                {
                    prho[facei] = tl.rhoPT(prho[facei], pp[facei], pT[facei]);
                }
                else
                {
                    prho[facei] =
                        tv.rhoPT(prho[facei], pp[facei], pT[facei])*xv
                      + tl.rhoPT(prho[facei], pp[facei], pT[facei])*xl;
                }
            }
        }
    }
}


template<class Thermo>
Foam::scalar Foam::cavitatingFluidBlastThermo<Thermo>::cellpRhoT
(
    const label celli,
    const bool limit
) const
{
    const scalar xv = this->cellx(celli);
    const scalar xl = 1.0 - xv;
    const scalar rho = this->rho_[celli];
    const scalar e = this->e_[celli];
    const scalar T = this->T_[celli];
    if (xv < this->residualFac_)
    {
        return Thermo::thermoType1::p(rho, e, T, limit);
    }
    else if (xl < this->residualFac_)
    {
        return Thermo::thermoType2::p(rho, e, T, limit);
    }

    return pSat_->pSat(T);
}


template<class Thermo>
Foam::scalar Foam::cavitatingFluidBlastThermo<Thermo>::patchFacepRhoT
(
    const label patchi,
    const label facei,
    const bool limit
) const
{
    const scalar xv = this->patchFacex(patchi, facei);
    const scalar xl = 1.0 - xv;
    const scalar rho = this->rho_.boundaryField()[patchi][facei];
    const scalar e = this->e_.boundaryField()[patchi][facei];
    const scalar T = this->T_.boundaryField()[patchi][facei];
    if (xv < this->residualFac_)
    {
        return Thermo::thermoType1::p(rho, e, T, limit);
    }
    else if (xl < this->residualFac_)
    {
        return Thermo::thermoType2::p(rho, e, T, limit);
    }

    return pSat_->pSat(T);
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::cavitatingFluidBlastThermo<Thermo>::Gamma() const
{
    tmp<volScalarField> tG
    (
        volScalarField::New
        (
            IOobject::groupName("Gamma", this->phaseName()),
            this->rho_.mesh(),
            1.0
        )
    );
    volScalarField& G = tG.ref();

    const typename Thermo::thermoType1& tl(*this);
    const typename Thermo::thermoType2& tv(*this);

    const scalarField& rhoI = this->rho_.primitiveField();
    const scalarField& eI = this->he().primitiveField();
    const scalarField& TI = this->T_.primitiveField();
    forAll(this->rho_, celli)
    {
        const scalar xv = this->cellx(celli);
        const scalar xl = 1.0 - xv;

        if (xl < this->residualFac_)
        {
            G[celli] = tv.Gamma(rhoI[celli], eI[celli], TI[celli]);
        }
        else if (xv < this->residualFac_)
        {
            G[celli] = tl.Gamma(rhoI[celli], eI[celli], TI[celli]);
        }
        else
        {
            G[celli] =
                tv.Gamma(rhoI[celli], eI[celli], TI[celli])*xv
              + tl.Gamma(rhoI[celli], eI[celli], TI[celli])*xl;
        }
    }

    volScalarField::Boundary& bG = G.boundaryFieldRef();

    forAll(bG, patchi)
    {
        scalarField& pG = bG[patchi];
        const scalarField& prho = this->rho_.boundaryField()[patchi];
        const scalarField& pe = this->he().boundaryField()[patchi];
        const scalarField& pT = this->T_.boundaryField()[patchi];
        const fvPatchScalarField& px = x_.boundaryField()[patchi];

        forAll(pG, facei)
        {
            const scalar xv = px[facei];
            const scalar xl = 1.0 - xv;
            if (xl < this->residualFac_)
            {
                pG[facei] = tv.Gamma(prho[facei], pe[facei], pT[facei]);
            }
            else if (xv < this->residualFac_)
            {
                pG[facei] = tl.Gamma(prho[facei], pe[facei], pT[facei]);
            }
            else
            {
                pG[facei] =
                    tv.Gamma(prho[facei], pe[facei], pT[facei])*xv
                  + tl.Gamma(prho[facei], pe[facei], pT[facei])*xl;
            }
        }
    }
    return tG;
}


template<class Thermo>
Foam::scalar
Foam::cavitatingFluidBlastThermo<Thermo>::cellGamma(const label celli) const
{
    const scalar x = this->cellx(celli);
    const scalar rho = this->rho_[celli];
    const scalar e = this->e_[celli];
    const scalar T = this->T_[celli];
    if (x < this->residualFac_)
    {
        return Thermo::thermoType1::Gamma(rho, e, T);
    }
    else if ((1.0 - x) < this->residualFac_)
    {
        return Thermo::thermoType2::Gamma(rho, e, T);
    }

    return
        Thermo::thermoType2::Gamma(rho, e, T)*x
      + Thermo::thermoType1::Gamma(rho, e, T)*(1.0 - x);
}


template<class Thermo>
Foam::scalar
Foam::cavitatingFluidBlastThermo<Thermo>::patchFaceGamma
(
    const label patchi,
    const label facei
) const
{
    const scalar x = this->patchFacex(patchi, facei);
    const scalar rho = this->rho_.boundaryField()[patchi][facei];
    const scalar e = this->e_.boundaryField()[patchi][facei];
    const scalar T = this->T_.boundaryField()[patchi][facei];
    if (x < this->residualFac_)
    {
        return Thermo::thermoType1::Gamma(rho, e, T);
    }
    else if ((1.0 - x) < this->residualFac_)
    {
        return Thermo::thermoType2::Gamma(rho, e, T);
    }

    return
        Thermo::thermoType2::Gamma(rho, e, T)*x
      + Thermo::thermoType1::Gamma(rho, e, T)*(1.0 - x);
}


template<class Thermo>
Foam::scalar
Foam::cavitatingFluidBlastThermo<Thermo>::celldpdRho(const label celli) const
{
    const scalar x = this->cellx(celli);
    const scalar rho = this->rho_[celli];
    const scalar e = this->e_[celli];
    const scalar T = this->T_[celli];
    if (x < this->residualFac_)
    {
        return Thermo::thermoType1::dpdRho(rho, e, T);
    }
    else if ((1.0 - x) < this->residualFac_)
    {
        return Thermo::thermoType2::dpdRho(rho, e, T);
    }

    return 0.0;
}


template<class Thermo>
Foam::scalar
Foam::cavitatingFluidBlastThermo<Thermo>::celldpde(const label celli) const
{
    const scalar x = this->cellx(celli);
    const scalar rho = this->rho_[celli];
    const scalar e = this->e_[celli];
    const scalar T = this->T_[celli];
    if (x < this->residualFac_)
    {
        return Thermo::thermoType1::dpde(rho, e, T);
    }
    else if ((1.0 - x) < this->residualFac_)
    {
        return Thermo::thermoType2::dpde(rho, e, T);
    }

    return 0.0;
}


template<class Thermo>
Foam::scalar
Foam::cavitatingFluidBlastThermo<Thermo>::celldpdT(const label celli) const
{
    const scalar x = this->cellx(celli);
    const scalar rho = this->rho_[celli];
    const scalar e = this->e_[celli];
    const scalar T = this->T_[celli];
    if (x < this->residualFac_)
    {
        return Thermo::thermoType1::dpdT(rho, e, T);
    }
    else if ((1.0 - x) < this->residualFac_)
    {
        return Thermo::thermoType2::dpdT(rho, e, T);
    }

    return pSat_->derivative(T);
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::cavitatingFluidBlastThermo<Thermo>::calce(const volScalarField& p) const
{
    tmp<volScalarField> teInit
    (
        volScalarField::New
        (
            IOobject::groupName("eInit", this->phaseName()),
            this->rho_.mesh(),
            dimensionedScalar(dimEnergy/dimMass, 0.0)
        )
    );
    volScalarField& eInit = teInit.ref();

    const typename Thermo::thermoType1& tl(*this);
    const typename Thermo::thermoType2& tv(*this);

    const scalarField& pI = p.primitiveField();
    const scalarField& rhoI = this->rho_.primitiveField();
    const scalarField& eI = this->he().primitiveField();
    const scalarField& TI = this->T_.primitiveField();
    forAll(this->rho_, celli)
    {
        const scalar xv = this->cellx(celli);
        const scalar xl = 1.0 - xv;

        if (xl < this->residualFac_)
        {
            eInit[celli] = tv.initializeEnergy
            (
                pI[celli],
                rhoI[celli],
                eI[celli],
                TI[celli]
            );
        }
        else if (xv < this->residualFac_)
        {
            eInit[celli] = tl.initializeEnergy
            (
                pI[celli],
                rhoI[celli],
                eI[celli],
                TI[celli]
            );
        }
        else
        {
            const scalar rhoSatv = rhoSatv_.lookup(TI[celli]);
            const scalar rhoSatl = rhoSatl_.lookup(TI[celli]);

            const scalar fv = rhoSatv/rhoI[celli];
            const scalar fl = rhoSatl/rhoI[celli];
            eInit[celli] =
                (
                    tv.initializeEnergy
                    (
                        pI[celli],
                        rhoI[celli],
                        eI[celli],
                        TI[celli]
                    )
                )*xv*fv
              + tl.initializeEnergy
                (
                    pI[celli],
                    rhoI[celli],
                    eI[celli],
                    this->T_[celli]
                )*xl*fl;
        }
    }

    volScalarField::Boundary& beInit = eInit.boundaryFieldRef();

    forAll(beInit, patchi)
    {
        scalarField& peInit = beInit[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& prho = this->rho_.boundaryField()[patchi];
        const scalarField& pe = this->he().boundaryField()[patchi];
        const scalarField& pT = this->T_.boundaryField()[patchi];
        const fvPatchScalarField& px = x_.boundaryField()[patchi];

        forAll(peInit, facei)
        {
            const scalar xv = px[facei];
            const scalar xl = 1.0 - xv;
            if (xl < this->residualFac_)
            {
                peInit[facei] = tv.initializeEnergy
                (
                    pp[facei],
                    prho[facei],
                    pe[facei],
                    pT[facei]
                );
            }
            else if (xv < this->residualFac_)
            {
                peInit[facei] = tl.initializeEnergy
                (
                    pp[facei],
                    prho[facei],
                    pe[facei],
                    pT[facei]
                );
            }
            else
            {
                const scalar rhoSatv = rhoSatv_.lookup(pT[facei]);
                const scalar rhoSatl = rhoSatl_.lookup(pT[facei]);

                const scalar fv = rhoSatv/prho[facei];
                const scalar fl = rhoSatl/prho[facei];
                peInit[facei] =
                    (
                        tv.initializeEnergy
                        (
                            pp[facei],
                            prho[facei],
                            pe[facei],
                            pT[facei]
                        )
                    )*xv*fv
                  + tl.initializeEnergy
                    (
                        pp[facei],
                        prho[facei],
                        pe[facei],
                        pT[facei]
                    )*xl*fl;
            }
        }
    }
    return teInit;
}


template<class Thermo>
Foam::scalar
Foam::cavitatingFluidBlastThermo<Thermo>::cellHE
(
    const scalar T,
    const label celli
) const
{
    const typename Thermo::thermoType1& tl(*this);
    const typename Thermo::thermoType2& tv(*this);

    const scalar rhoi = this->rho_[celli];
    const scalar ei = this->e_[celli];
    const scalar xv = this->cellx(celli);
    const scalar xl = 1.0 - xv;

    if (xl < small)
    {
        return tv.Es(rhoi, ei, T);
    }
    else if (xv < small)
    {
        return tl.Es(rhoi, ei, T);
    }
    else if (rhoi > small)
    {
        const scalar rhoSatv = rhoSatv_.lookup(T);
        const scalar rhoSatl = rhoSatl_.lookup(T);

        const scalar fv = rhoSatv/rhoi;
        const scalar fl = rhoSatl/rhoi;

        return
            tv.Es(rhoi, ei, T)*xv*fv
          + tl.Es(rhoi, ei, T)*xl*fl;
    }
    return 0.0;
}


template<class Thermo>
Foam::scalar
Foam::cavitatingFluidBlastThermo<Thermo>::patchFaceHE
(
    const scalar T,
    const label patchi,
    const label facei
) const
{
    const typename Thermo::thermoType1& tl(*this);
    const typename Thermo::thermoType2& tv(*this);

    const scalar rhoi = this->rho_.boundaryField()[patchi][facei];
    const scalar ei = this->e_.boundaryField()[patchi][facei];
    const scalar xv = this->patchFacex(patchi, facei);
    const scalar xl = 1.0 - xv;

    if (xl < small)
    {
        return tv.Es(rhoi, ei, T);
    }
    else if (xv < small)
    {
        return tl.Es(rhoi, ei, T);
    }
    else if (rhoi > small)
    {
        const scalar rhoSatv = rhoSatv_.lookup(T);
        const scalar rhoSatl = rhoSatl_.lookup(T);

        const scalar fv = rhoSatv/rhoi;
        const scalar fl = rhoSatl/rhoi;

        return
            tv.Es(rhoi, ei, T)*xv*fv
          + tl.Es(rhoi, ei, T)*xl*fl;
    }
    return 0.0;
}


template<class Thermo>
Foam::scalar
Foam::cavitatingFluidBlastThermo<Thermo>::cellTHE
(
    const scalar he,
    const scalar T0,
    const label celli
) const
{
    const typename Thermo::thermoType1& tl(*this);
    const typename Thermo::thermoType2& tv(*this);

    const scalar rhoi = this->rho_[celli];
    const scalar xv = this->cellx(celli);
    const scalar xl = 1.0 - xv;

    if (xl < small)
    {
        return tv.TRhoE(T0, rhoi, he);
    }
    else if (xv < small)
    {
        return tl.TRhoE(T0, rhoi, he);
    }
    else
    {
        return
            tv.TRhoE(T0, rhoi, he)*xv
          + tl.TRhoE(T0, rhoi, he)*xl;
    }
}


template<class Thermo>
Foam::scalar
Foam::cavitatingFluidBlastThermo<Thermo>::calcCelle
(
    const scalar p,
    const label celli
) const
{
    const scalar xv = this->cellx(celli);
    const scalar xl = 1.0 - xv;

    const scalar rhoi = this->rho_[celli];
    const scalar e = this->e_[celli];
    const scalar T = this->T_[celli];
    if (xv < this->residualFac_)
    {
        return Thermo::thermoType1::initializeEnergy(p, rhoi, e, T);
    }
    else if (xl < this->residualFac_)
    {
        return Thermo::thermoType2::initializeEnergy(p, rhoi, e, T);
    }

    const scalar rhoSatv = rhoSatv_.lookup(T);
    const scalar rhoSatl = rhoSatl_.lookup(T);

    const scalar fv = rhoSatv/rhoi;
    const scalar fl = rhoSatl/rhoi;

    return
        Thermo::thermoType2::initializeEnergy(p, rhoi, e, T)*xv*fv
      + Thermo::thermoType1::initializeEnergy(p, rhoi, e, T)*xl*fl;
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::cavitatingFluidBlastThermo<Thermo>::ESource() const
{
    return tmp<volScalarField>
    (
        volScalarField::New
        (
            "ESource",
            dxdt_()*hv_*this->rho_
            // this->rho_.mesh(),
            // dimensionedScalar("0", dimEnergy/dimTime/dimVolume, 0.0)
        )
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::cavitatingFluidBlastThermo<Thermo>::initESource() const
{
    return volScalarField::New
    (
        "initESource",
        // x_*hv_
        this->rho_.mesh(),
        dimensionedScalar("0", dimEnergy/dimMass, 0.0)
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::cavitatingFluidBlastThermo<Thermo>::calcp() const
{
   tmp<volScalarField> tp
    (
        volScalarField::New
        (
            IOobject::groupName("p", this->phaseName()),
            this->rho_.mesh(),
            dimensionedScalar(dimPressure, 0.0)
        )
    );
    volScalarField& p = tp.ref();

    const typename Thermo::thermoType1& tl(*this);
    const typename Thermo::thermoType2& tv(*this);
    forAll(this->rho_, celli)
    {
        const scalar rhoi = this->rho_[celli];
        const scalar ei = this->e_[celli];
        const scalar Ti = this->T_[celli];

        const scalar xv = this->cellx(celli);
        const scalar xl = 1.0 - xv;

        if (xl < this->residualFac_)
        {
            p[celli] = tv.pRhoT(rhoi, ei, Ti);
        }
        else if (xv < this->residualFac_)
        {
            p[celli] = tl.pRhoT(rhoi, ei, Ti);
        }
        else
        {
            p[celli] = pSat_->pSat(Ti);
        }
    }

    volScalarField::Boundary& bp = p.boundaryFieldRef();
    forAll(bp, patchi)
    {
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        const fvPatchScalarField& phe = this->e_.boundaryField()[patchi];
        const fvPatchScalarField& px = x_.boundaryField()[patchi];

        fvPatchScalarField& pp = bp[patchi];

        forAll(pp, facei)
        {
            const scalar rhoi = prho[facei];
            const scalar ei = phe[facei];
            const scalar Ti = pT[facei];

            const scalar xv = px[facei];
            const scalar xl = 1.0 - xv;


            if (xl < this->residualFac_)
            {
                pp[facei] = tv.pRhoT(rhoi, ei, Ti);
            }
            else if (xv < this->residualFac_)
            {
                pp[facei] = tl.pRhoT(rhoi, ei, Ti);
            }
            else
            {
                pp[facei] = pSat_->pSat(Ti);
            }
        }
    }
    return tp;
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::cavitatingFluidBlastThermo<Thermo>::calcSpeedOfSound() const
{
    tmp<volScalarField> tc
    (
        volScalarField::New
        (
            "cSqr",
            this->rho_.mesh(),
            dimensionedScalar(dimVelocity, 0.0)
        )
    );
    volScalarField& c = tc.ref();

    const typename Thermo::thermoType1& tl(*this);
    const typename Thermo::thermoType2& tv(*this);
    forAll(this->rho_, celli)
    {
        const scalar rhoi = this->rho_[celli];
        const scalar ei = this->e_[celli];
        const scalar Ti = this->T_[celli];
        const scalar pi = this->p_[celli];

        const scalar xv = this->cellx(celli);
        const scalar xl = 1.0 - xv;

        if (xl < this->residualFac_)
        {
            c[celli] = sqrt(max(tv.cSqr(pi, rhoi, ei, Ti), small));
        }
        else if (xv < this->residualFac_)
        {
            c[celli] = sqrt(max(tl.cSqr(pi, rhoi, ei, Ti), small));
        }
        else
        {
            const scalar rhoSatv = rhoSatv_.lookup(Ti);
            const scalar rhoSatl = rhoSatl_.lookup(Ti);

            c[celli] =
                sqrt
                (
                    1.0
                   /(
                        xv/(rhoSatv*max(tv.cSqr(pi, rhoi, ei, Ti), small))
                      + xl/(rhoSatl*max(tl.cSqr(pi, rhoi, ei, Ti), small))
                    )/rhoi
                );
        }
    }

    volScalarField::Boundary& bc = c.boundaryFieldRef();
    forAll(bc, patchi)
    {
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        const fvPatchScalarField& phe = this->e_.boundaryField()[patchi];
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& px = x_.boundaryField()[patchi];

        fvPatchScalarField& pc = bc[patchi];

        forAll(pc, facei)
        {
            const scalar rhoi = prho[facei];
            const scalar ei = phe[facei];
            const scalar Ti = pT[facei];
            const scalar pi = pp[facei];

            const scalar xv = px[facei];
            const scalar xl = 1.0 - xv;


            if (xl < this->residualFac_)
            {
                pc[facei] = sqrt(max(tv.cSqr(pi, rhoi, ei, Ti), small));
            }
            else if (xv < this->residualFac_)
            {
                pc[facei] = sqrt(max(tl.cSqr(pi, rhoi, ei, Ti), small));
            }
            else
            {
                const scalar rhoSatv = rhoSatl_.lookup(Ti);
                const scalar rhoSatl = rhoSatv_.lookup(Ti);
                pc[facei] =
                    sqrt
                    (
                        1.0
                       /(
                            xv/(rhoSatv*max(tv.cSqr(pi, rhoi, ei, Ti), small))
                          + xl/(rhoSatl*max(tl.cSqr(pi, rhoi, ei, Ti), small))
                        )/rhoi
                    );
            }
        }
    }
    return tc;
}


// ************************************************************************* //
