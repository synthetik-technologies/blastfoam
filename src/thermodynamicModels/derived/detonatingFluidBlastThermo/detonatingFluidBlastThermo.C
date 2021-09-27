/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a detonating material
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

#include "detonatingFluidBlastThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
void Foam::detonatingFluidBlastThermo<Thermo>::calculate()
{
    const typename Thermo::thermoType1& t1(*this);
    const typename Thermo::thermoType2& t2(*this);
    forAll(this->rho_, celli)
    {
        const scalar x2 = this->cellx(celli);
        const scalar x1 = 1.0 - x2;
        const scalar rhoi(this->rho_[celli]);
        scalar& ei(this->he()[celli]);
        scalar& Ti(this->Thermo::baseThermo::T()[celli]);

        if (x2 < small)
        {
            Ti =
                t1.TRhoE(Ti, rhoi, ei);
            if (Ti < this->TLow_)
            {
                ei = t1.Es(rhoi, ei, this->TLow_);
                Ti = this->TLow_;
            }

            const scalar pi = t1.p(rhoi, ei, Ti);
            const scalar Cpi = t1.Cp(rhoi, ei, Ti);

            this->Thermo::baseThermo::p()[celli] = pi;
            this->CpRef()[celli] = Cpi;
            this->CvRef()[celli] = t1.Cv(rhoi, ei, Ti);
            this->muRef()[celli] = t1.mu(rhoi, ei, Ti);
            this->alpha()[celli] = t1.kappa(rhoi, ei, Ti)/Cpi;
            this->speedOfSound()[celli] =
                sqrt(max(t1.cSqr(pi, rhoi, ei, Ti), small));
        }
        else if (x1 < small)
        {
            Ti =
                t2.TRhoE(Ti, rhoi, ei);
            if (Ti < this->TLow_)
            {
                ei = t2.Es(rhoi, ei, this->TLow_);
                Ti = this->TLow_;
            }

            const scalar pi = t2.p(rhoi, ei, Ti);
            const scalar Cpi = t2.Cp(rhoi, ei, Ti);

            this->Thermo::baseThermo::p()[celli] = pi;
            this->CpRef()[celli] = Cpi;
            this->CvRef()[celli] = t2.Cv(rhoi, ei, Ti);
            this->muRef()[celli] = t2.mu(rhoi, ei, Ti);
            this->alpha()[celli] = t2.kappa(rhoi, ei, Ti)/Cpi;
            this->speedOfSound()[celli] =
                sqrt(max(t2.cSqr(pi, rhoi, ei, Ti), small));
        }
        else
        {
            Ti =
                t1.TRhoE(Ti, rhoi, ei)*x1
              + t2.TRhoE(Ti, rhoi, ei)*x2;
            if (Ti < this->TLow_)
            {
                ei =
                    t1.Es(rhoi, ei, this->TLow_)*x1
                  + t2.Es(rhoi, ei, this->TLow_)*x2;
                Ti = this->TLow_;
            }

            const scalar pi =
                t1.p(rhoi, ei, Ti)*x1
              + t2.p(rhoi, ei, Ti)*x2;

            this->Thermo::baseThermo::p()[celli] = pi;
            this->CpRef()[celli] =
                t1.Cp(rhoi, ei, Ti)*x1
              + t2.Cp(rhoi, ei, Ti)*x2;;
            this->CvRef()[celli] =
                t1.Cv(rhoi, ei, Ti)*x1
              + t2.Cv(rhoi, ei, Ti)*x2;
            this->muRef()[celli] =
                t1.mu(rhoi, ei, Ti)*x1
              + t2.mu(rhoi, ei, Ti)*x2;
            this->alpha()[celli] =
                t1.kappa(rhoi, ei, Ti)/t1.Cp(rhoi, ei, Ti)*x1
              + t2.kappa(rhoi, ei, Ti)/t2.Cp(rhoi, ei, Ti)*x2;
            this->speedOfSound()[celli] =
                sqrt
                (
                    max(t1.cSqr(pi, rhoi, ei, Ti), small)*x1
                  + max(t2.cSqr(pi, rhoi, ei, Ti), small)*x2
                );
        }
    }

    this->Thermo::baseThermo::T().correctBoundaryConditions();
    this->he().correctBoundaryConditions();
    this->Thermo::baseThermo::p().correctBoundaryConditions();

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pT =
            this->Thermo::baseThermo::T().boundaryField()[patchi];
        const fvPatchScalarField& phe = this->he().boundaryField()[patchi];
        const fvPatchScalarField& pp =
            this->Thermo::baseThermo::p().boundaryField()[patchi];
        const scalarField px(this->x(patchi));

        fvPatchScalarField& pCp = this->CpRef().boundaryFieldRef()[patchi];
        fvPatchScalarField& pCv = this->CvRef().boundaryFieldRef()[patchi];
        fvPatchScalarField& pmu = this->muRef().boundaryFieldRef()[patchi];
        fvPatchScalarField& palpha =
            this->alpha().boundaryFieldRef()[patchi];
        fvPatchScalarField& pc =
            this->speedOfSound().boundaryFieldRef()[patchi];

        forAll(pT, facei)
        {
            const scalar x2 = px[facei];
            const scalar x1 = 1.0 - x2;
            const scalar rhoi(prho[facei]);
            const scalar ei(phe[facei]);
            const scalar Ti(pT[facei]);
            const scalar pi(pp[facei]);

            if (x2 < small)
            {
                pCp[facei] = t1.Cp(rhoi, ei, Ti);
                pCv[facei] = t1.Cv(rhoi, ei, Ti);
                pmu[facei] = t1.mu(rhoi, ei, Ti);
                palpha[facei] = t1.kappa(rhoi, ei, Ti)/pCp[facei];
                pc[facei] =
                    sqrt(max(t1.cSqr(pi, rhoi, ei, Ti), small));
            }
            else if (x1 < small)
            {
                pCp[facei] = t2.Cp(rhoi, ei, Ti);
                pCv[facei] = t2.Cv(rhoi, ei, Ti);
                pmu[facei] = t2.mu(rhoi, ei, Ti);
                palpha[facei] = t2.kappa(rhoi, ei, Ti)/pCp[facei];
                pc[facei] =
                    sqrt(max(t2.cSqr(pi, rhoi, ei, Ti), small));
            }
            else
            {
                pCp[facei] =
                    t1.Cp(rhoi, ei, Ti)*x1
                  + t2.Cp(rhoi, ei, Ti)*x2;
                pCv[facei] =
                    t1.Cv(rhoi, ei, Ti)*x1
                  + t2.Cv(rhoi, ei, Ti)*x2;
                pmu[facei] =
                    t1.mu(rhoi, ei, Ti)*x1
                  + t2.mu(rhoi, ei, Ti)*x2;
                palpha[facei] =
                    t1.kappa(rhoi, ei, Ti)/pCp[facei]*x1
                  + t2.kappa(rhoi, ei, Ti)/pCp[facei]*x2;
                pc[facei] =
                    sqrt
                    (
                        max(t1.cSqr(pi, rhoi, ei, Ti), small)*x1
                      + max(t2.cSqr(pi, rhoi, ei, Ti), small)*x2
                    );
            }
        }
    }
}



template<class Thermo>
void Foam::detonatingFluidBlastThermo<Thermo>::calculate
(
    const volScalarField& alpha,
    const volScalarField& he,
    const volScalarField& T,
    volScalarField& alphaCp,
    volScalarField& alphaCv,
    volScalarField& alphaMu,
    volScalarField& alphaAlphah,
    volScalarField& pXiSum,
    volScalarField& XiSum,
    volScalarField& rhoXiSum
)
{
    const typename Thermo::thermoType1& t1(*this);
    const typename Thermo::thermoType2& t2(*this);

    forAll(alpha, celli)
    {
        const scalar x2 = this->cellx(celli);
        const scalar x1 = 1.0 - x2;
        const scalar alphai = alpha[celli];
        const scalar rhoi(this->rho_[celli]);
        const scalar ei(he[celli]);
        const scalar Ti(T[celli]);
        if (alphai > this->residualAlpha().value())
        {
            scalar Xii = alphai;

            if (x2 < small)
            {
                Xii /= t1.Gamma(rhoi, ei, Ti) - 1.0;

                alphaCp[celli] += t1.Cp(rhoi, ei, Ti)*alphai;
                alphaCv[celli] += t1.Cv(rhoi, ei, Ti)*alphai;
                alphaMu[celli] += t1.mu(rhoi, ei, Ti)*alphai;
                alphaAlphah[celli] +=
                    t1.kappa(rhoi, ei, Ti)/t1.Cp(rhoi, ei, Ti)*alphai;
                pXiSum[celli] += t1.p(rhoi, ei, Ti)*Xii;
                XiSum[celli] += Xii;
                rhoXiSum[celli] += rhoi*Xii;
            }
            else if (x1 < small)
            {
                Xii /= t2.Gamma(rhoi, ei, Ti) - 1.0;

                alphaCp[celli] += t2.Cp(rhoi, ei, Ti)*alphai;
                alphaCv[celli] += t2.Cv(rhoi, ei, Ti)*alphai;
                alphaMu[celli] += t2.mu(rhoi, ei, Ti)*alphai;
                alphaAlphah[celli] +=
                    t2.kappa(rhoi, ei, Ti)/t2.Cp(rhoi, ei, Ti)*alphai;
                pXiSum[celli] += t2.p(rhoi, ei, Ti)*Xii;
                XiSum[celli] += Xii;
                rhoXiSum[celli] += rhoi*Xii;
            }
            else
            {
                Xii /=
                    x1*t1.Gamma(rhoi, ei, Ti)
                  + x2*t2.Gamma(rhoi, ei, Ti)
                  - 1.0;
                alphaCp[celli] +=
                    (
                        t1.Cp(rhoi, ei, Ti)*x1
                      + t2.Cp(rhoi, ei, Ti)*x2
                    )*alphai;
                alphaCv[celli] +=
                    (
                        t1.Cv(rhoi, ei, Ti)*x1
                      + t2.Cv(rhoi, ei, Ti)*x2
                    )*alphai;
                alphaMu[celli] +=
                    (
                        t1.mu(rhoi, ei, Ti)*x1
                      + t2.mu(rhoi, ei, Ti)*x2
                    )*alphai;
                alphaAlphah[celli] +=
                    (
                        t1.kappa(rhoi, ei, Ti)/t1.Cp(rhoi, ei, Ti)*x1
                      + t2.kappa(rhoi, ei, Ti)/t2.Cp(rhoi, ei, Ti)*x2
                    )*alphai;
                pXiSum[celli] +=
                    (
                        t1.p(rhoi, ei, Ti)*x1
                      + t2.p(rhoi, ei, Ti)*x2
                    )*Xii;
                XiSum[celli] += Xii;
                rhoXiSum[celli] += rhoi*Xii;
            }
        }
    }

    forAll(alpha.boundaryField(), patchi)
    {
        const fvPatchScalarField& palpha = alpha.boundaryField()[patchi];
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pT = T.boundaryField()[patchi];
        const fvPatchScalarField& phe = he.boundaryField()[patchi];
        const scalarField px(this->x(patchi));

        fvPatchScalarField& palphaCp = alphaCp.boundaryFieldRef()[patchi];
        fvPatchScalarField& palphaCv = alphaCv.boundaryFieldRef()[patchi];
        fvPatchScalarField& palphaMu = alphaMu.boundaryFieldRef()[patchi];
        fvPatchScalarField& palphaAlphah =
            alphaAlphah.boundaryFieldRef()[patchi];
        fvPatchScalarField& ppXiSum = pXiSum.boundaryFieldRef()[patchi];
        fvPatchScalarField& pxiSum = XiSum.boundaryFieldRef()[patchi];
        fvPatchScalarField& prhoXiSum = rhoXiSum.boundaryFieldRef()[patchi];

        forAll(palpha, facei)
        {

            const scalar x2 = px[facei];
            const scalar x1 = 1.0 - x2;
            const scalar alphai = palpha[facei];
            const scalar& rhoi(prho[facei]);
            const scalar& ei(phe[facei]);
            const scalar& Ti(pT[facei]);
            scalar Xii = alphai;

            if (x2 < small)
            {
                Xii /= t1.Gamma(rhoi, ei, Ti) - 1.0;

                ppXiSum[facei] += t1.p(rhoi, ei, Ti)*Xii;
                palphaCp[facei] += t1.Cp(rhoi, ei, Ti)*alphai;
                palphaCv[facei] += t1.Cv(rhoi, ei, Ti)*alphai;
                palphaMu[facei] += t1.mu(rhoi, ei, Ti)*alphai;
                palphaAlphah[facei] +=
                    t1.kappa(rhoi, ei, Ti)/t1.Cp(rhoi, ei, Ti)*alphai;
            }
            else if (x1 < small)
            {
                Xii /= t2.Gamma(rhoi, ei, Ti) - 1.0;

                ppXiSum[facei] += t2.p(rhoi, ei, Ti);
                palphaCp[facei] += t2.Cp(rhoi, ei, Ti);
                palphaCv[facei] += t2.Cv(rhoi, ei, Ti);
                palphaMu[facei] += t2.mu(rhoi, ei, Ti);
                palphaAlphah[facei] +=
                    t2.kappa(rhoi, ei, Ti)/t2.Cp(rhoi, ei, Ti);
            }
            else
            {
                Xii /=
                    t1.Gamma(rhoi, ei, Ti)*x1
                  + t2.Gamma(rhoi, ei, Ti)*x2
                  - 1.0;
                ppXiSum[facei] +=
                    t1.p(rhoi, ei, Ti)*x1
                  + t2.p(rhoi, ei, Ti)*x2;
                palphaCp[facei] +=
                    t1.Cp(rhoi, ei, Ti)*x1
                  + t2.Cp(rhoi, ei, Ti)*x2;
                palphaCv[facei] +=
                    t1.Cv(rhoi, ei, Ti)*x1
                  + t2.Cv(rhoi, ei, Ti)*x2;
                palphaMu[facei] +=
                    t1.mu(rhoi, ei, Ti)*x1
                  + t2.mu(rhoi, ei, Ti)*x2;
                palphaAlphah[facei] +=
                    t1.kappa(rhoi, ei, Ti)/t1.Cp(rhoi, ei, Ti)*x1
                  + t2.kappa(rhoi, ei, Ti)/t2.Cp(rhoi, ei, Ti)*x2;
            }
            pxiSum[facei] += Xii;
            prhoXiSum[facei] += rhoi*Xii;
        }
    }
}


template<class Thermo>
void Foam::detonatingFluidBlastThermo<Thermo>::calculateSpeedOfSound
(
    const volScalarField& alpha,
    volScalarField& cSqrRhoXiSum
)
{
    const typename Thermo::thermoType1& t1(*this);
    const typename Thermo::thermoType2& t2(*this);

    forAll(this->rho_, celli)
    {
        if (alpha[celli] > this->residualAlpha().value())
        {
            const scalar x2 = this->cellx(celli);
            const scalar x1 = 1.0 - x2;

            if (x2 < small)
            {
                cSqrRhoXiSum[celli] +=
                    t1.cSqr
                    (
                        this->p_[celli],
                        this->rho_[celli],
                        this->e_[celli],
                        this->T_[celli]
                    )*this->rho_[celli]*alpha[celli]
                   /(
                        t1.Gamma
                        (
                            this->rho_[celli],
                            this->e_[celli],
                            this->T_[celli]
                        )
                      - 1.0
                    );
            }
            else if (x1 < small)
            {
                cSqrRhoXiSum[celli] +=
                    t2.cSqr
                    (
                        this->p_[celli],
                        this->rho_[celli],
                        this->e_[celli],
                        this->T_[celli]
                    )*this->rho_[celli]*alpha[celli]
                   /(
                        t2.Gamma
                        (
                            this->rho_[celli],
                            this->e_[celli],
                            this->T_[celli]
                        )
                      - 1.0
                    );
            }
            else
            {
                cSqrRhoXiSum[celli] +=
                    (
                        t1.cSqr
                        (
                            this->p_[celli],
                            this->rho_[celli],
                            this->e_[celli],
                            this->T_[celli]
                        )*x1
                      + t2.cSqr
                        (
                            this->p_[celli],
                            this->rho_[celli],
                            this->e_[celli],
                            this->T_[celli]
                        )*x2
                    )*this->rho_[celli]*alpha[celli]
                   /(
                        t1.Gamma
                        (
                            this->rho_[celli],
                            this->e_[celli],
                            this->T_[celli]
                        )*x1
                      + t2.Gamma
                        (
                            this->rho_[celli],
                            this->e_[celli],
                            this->T_[celli]
                        )  *x2
                      - 1.0
                    );
            }
        }
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& palpha = alpha.boundaryField()[patchi];
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        const fvPatchScalarField& phe = this->e_.boundaryField()[patchi];
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const scalarField px(this->x(patchi));
        fvPatchScalarField& pcSqr = cSqrRhoXiSum.boundaryFieldRef()[patchi];

        forAll(pT, facei)
        {
            const scalar x2 = px[facei];
            const scalar x1 = 1.0 - x2;
            if (x2 < small)
            {
                pcSqr[facei] +=
                    t1.cSqr
                    (
                        pp[facei],
                        prho[facei],
                        phe[facei],
                        pT[facei]
                    )*prho[facei]*palpha[facei]
                   /(
                        t1.Gamma
                        (
                            prho[facei],
                            phe[facei],
                            pT[facei]
                        )
                      - 1.0
                    );
            }
            else if (x1 < small)
            {
                pcSqr[facei] +=
                    t2.cSqr
                    (
                        pp[facei],
                        prho[facei],
                        phe[facei],
                        pT[facei]
                    )*prho[facei]*palpha[facei]
                   /(
                        t2.Gamma
                        (
                            prho[facei],
                            phe[facei],
                            pT[facei]
                        )
                      - 1.0
                    );
            }
            else
            {
                pcSqr[facei] +=
                    (
                        t1.cSqr
                        (
                            pp[facei],
                            prho[facei],
                            phe[facei],
                            pT[facei]
                        )*x1
                      + t2.cSqr
                        (
                            pp[facei],
                            prho[facei],
                            phe[facei],
                            pT[facei]
                        )*x2
                    )*prho[facei]*palpha[facei]
                   /(
                        t1.Gamma
                        (
                            prho[facei],
                            phe[facei],
                            pT[facei]
                        )*x1
                      + t2.Gamma
                        (
                            prho[facei],
                            phe[facei],
                            pT[facei]
                        )  *x2
                      - 1.0
                    );
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::detonatingFluidBlastThermo<Thermo>::detonatingFluidBlastThermo
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
        dict.subDict("reactants"),
        dict.subDict("products"),
        phaseName,
        masterName
    ),
    activation_
    (
        activationModel::New
        (
            mesh,
            dict,
            phaseName
        )
    ),
    afterburn_
    (
        afterburnModel::New
        (
            mesh,
            dict,
            phaseName
        )
    )
{
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
void Foam::detonatingFluidBlastThermo<Thermo>::initializeModels()
{
    activation_->initializeModels();
    afterburn_->initializeModels();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::detonatingFluidBlastThermo<Thermo>::~detonatingFluidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::detonatingFluidBlastThermo<Thermo>::correct()
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
void Foam::detonatingFluidBlastThermo<Thermo>::update()
{
    activation_->update();
    afterburn_->update();
}


template<class Thermo>
void Foam::detonatingFluidBlastThermo<Thermo>::solve()
{
    activation_->solve();
    afterburn_->solve();
}


template<class Thermo>
void Foam::detonatingFluidBlastThermo<Thermo>::postUpdate()
{
    activation_->postUpdate();
    afterburn_->postUpdate();
}


template<class Thermo>
void Foam::detonatingFluidBlastThermo<Thermo>::updateRho(const volScalarField& p)
{
    volScalarField rhoNew
    (
        Thermo::blendedVolScalarFieldProperty
        (
            "rho",
            dimDensity,
            &Thermo::thermoType1::rhoPT,
            &Thermo::thermoType2::rhoPT,
            p,
            this->T_
        )
    );
    this->rho_ = rhoNew;
    forAll(this->rho_.boundaryField(), patchi)
    {
        forAll(this->rho_.boundaryField()[patchi], facei)
        {
            this->rho_.boundaryFieldRef()[patchi][facei] =
                rhoNew.boundaryField()[patchi][facei];
        }
    }
}


template<class Thermo>
Foam::scalar
Foam::detonatingFluidBlastThermo<Thermo>::cellpRhoT(const label celli) const
{
    const scalar& x = this->cellx(celli);
    const scalar rho = this->rho_[celli];
    const scalar e = this->e_[celli];
    const scalar T = this->T_[celli];
    if (x < small)
    {
        return Thermo::thermoType1::p(rho, e, T);
    }
    else if ((1.0 - x) < small)
    {
        return Thermo::thermoType2::p(rho, e, T);
    }

    return
        Thermo::thermoType2::p(rho, e, T)*x
      + Thermo::thermoType1::p(rho, e, T)*(1.0 - x);
}


template<class Thermo>
Foam::scalar
Foam::detonatingFluidBlastThermo<Thermo>::cellGamma(const label celli) const
{
    const scalar& x = this->cellx(celli);
    const scalar rho = this->rho_[celli];
    const scalar e = this->e_[celli];
    const scalar T = this->T_[celli];
    if (x < small)
    {
        return Thermo::thermoType1::Gamma(rho, e, T);
    }
    else if ((1.0 - x) < small)
    {
        return Thermo::thermoType2::Gamma(rho, e, T);
    }

    return
        Thermo::thermoType2::Gamma(rho, e, T)*x
      + Thermo::thermoType1::Gamma(rho, e, T)*(1.0 - x);
}


template<class Thermo>
Foam::scalar
Foam::detonatingFluidBlastThermo<Thermo>::celldpdRho(const label celli) const
{
    const scalar& x = this->cellx(celli);
    const scalar rho = this->rho_[celli];
    const scalar e = this->e_[celli];
    const scalar T = this->T_[celli];
    if (x < small)
    {
        return Thermo::thermoType1::dpdRho(rho, e, T);
    }
    else if ((1.0 - x) < small)
    {
        return Thermo::thermoType2::dpdRho(rho, e, T);
    }

    return
        Thermo::thermoType2::dpdRho(rho, e, T)*x
      + Thermo::thermoType1::dpdRho(rho, e, T)*(1.0 - x);
}


template<class Thermo>
Foam::scalar
Foam::detonatingFluidBlastThermo<Thermo>::celldpde(const label celli) const
{
    const scalar& x = this->cellx(celli);
    const scalar rho = this->rho_[celli];
    const scalar e = this->e_[celli];
    const scalar T = this->T_[celli];
    if (x < small)
    {
        return Thermo::thermoType1::dpde(rho, e, T);
    }
    else if ((1.0 - x) < small)
    {
        return Thermo::thermoType2::dpde(rho, e, T);
    }

    return
        Thermo::thermoType2::dpde(rho, e, T)*x
      + Thermo::thermoType1::dpde(rho, e, T)*(1.0 - x);
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidBlastThermo<Thermo>::calce(const volScalarField& p) const
{
    tmp<volScalarField> eInit
    (
        Thermo::volScalarFieldProperty
        (
            "e",
            dimEnergy/dimMass,
            &Thermo::thermoType1::initializeEnergy,
            p,
            this->rho_,
            this->e_,
            this->T_
        )
    );

    //- Add detonation energy to initially reacted material
    if (this->rho_.time().timeIndex() == 0)
    {
        eInit.ref() += activation_->e0()*activation_->lambda();
    }

    return eInit;
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidBlastThermo<Thermo>::ESource() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "ESource",
                this->rho_.mesh().time().timeName(),
                this->rho_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            (activation_->ESource() + afterburn_->ESource())*this->rho_
        )
    );
}



// ************************************************************************* //
