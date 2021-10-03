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
        scalar& ei(this->heRef()[celli]);
        scalar& Ti(this->TRef()[celli]);

        if (x2 < this->residualActivation_)
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

            this->pRef()[celli] = pi;
            this->CpRef()[celli] = Cpi;
            this->CvRef()[celli] = t1.Cv(rhoi, ei, Ti);
            this->muRef()[celli] = t1.mu(rhoi, ei, Ti);
            this->alphaRef()[celli] = t1.kappa(rhoi, ei, Ti)/Cpi;
            this->speedOfSoundRef()[celli] =
                sqrt(max(t1.cSqr(pi, rhoi, ei, Ti), small));
        }
        else if (x1 < this->residualActivation_)
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

            this->pRef()[celli] = pi;
            this->CpRef()[celli] = Cpi;
            this->CvRef()[celli] = t2.Cv(rhoi, ei, Ti);
            this->muRef()[celli] = t2.mu(rhoi, ei, Ti);
            this->alphaRef()[celli] = t2.kappa(rhoi, ei, Ti)/Cpi;
            this->speedOfSoundRef()[celli] =
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

            this->pRef()[celli] = pi;
            this->CpRef()[celli] =
                t1.Cp(rhoi, ei, Ti)*x1
              + t2.Cp(rhoi, ei, Ti)*x2;;
            this->CvRef()[celli] =
                t1.Cv(rhoi, ei, Ti)*x1
              + t2.Cv(rhoi, ei, Ti)*x2;
            this->muRef()[celli] =
                t1.mu(rhoi, ei, Ti)*x1
              + t2.mu(rhoi, ei, Ti)*x2;
            this->alphaRef()[celli] =
                t1.kappa(rhoi, ei, Ti)/t1.Cp(rhoi, ei, Ti)*x1
              + t2.kappa(rhoi, ei, Ti)/t2.Cp(rhoi, ei, Ti)*x2;
            this->speedOfSoundRef()[celli] =
                sqrt
                (
                    max(t1.cSqr(pi, rhoi, ei, Ti), small)*x1
                  + max(t2.cSqr(pi, rhoi, ei, Ti), small)*x2
                );
        }
    }

    this->TRef().correctBoundaryConditions();
    this->heRef().correctBoundaryConditions();
    this->pRef().correctBoundaryConditions();

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pT =
            this->TRef().boundaryField()[patchi];
        const fvPatchScalarField& phe =
            this->heRef().boundaryField()[patchi];
        const fvPatchScalarField& pp =
            this->pRef().boundaryField()[patchi];
        const scalarField px(this->x(patchi));

        fvPatchScalarField& pCp = this->CpRef().boundaryFieldRef()[patchi];
        fvPatchScalarField& pCv = this->CvRef().boundaryFieldRef()[patchi];
        fvPatchScalarField& pmu = this->muRef().boundaryFieldRef()[patchi];
        fvPatchScalarField& palpha =
            this->alphaRef().boundaryFieldRef()[patchi];
        fvPatchScalarField& pc =
            this->speedOfSoundRef().boundaryFieldRef()[patchi];

        forAll(pT, facei)
        {
            const scalar x2 = px[facei];
            const scalar x1 = 1.0 - x2;
            const scalar rhoi(prho[facei]);
            const scalar ei(phe[facei]);
            const scalar Ti(pT[facei]);
            const scalar pi(pp[facei]);

            if (x2 < this->residualActivation_)
            {
                pCp[facei] = t1.Cp(rhoi, ei, Ti);
                pCv[facei] = t1.Cv(rhoi, ei, Ti);
                pmu[facei] = t1.mu(rhoi, ei, Ti);
                palpha[facei] = t1.kappa(rhoi, ei, Ti)/pCp[facei];
                pc[facei] =
                    sqrt(max(t1.cSqr(pi, rhoi, ei, Ti), small));
            }
            else if (x1 < this->residualActivation_)
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
    volScalarField& XiSum
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
        if (alphai > this->residualAlpha_.value())
        {
            scalar Gamma = alphai;
            scalar pi;

            if (x2 < this->residualActivation_)
            {
                alphaCp[celli] += t1.Cp(rhoi, ei, Ti)*alphai;
                alphaCv[celli] += t1.Cv(rhoi, ei, Ti)*alphai;
                alphaMu[celli] += t1.mu(rhoi, ei, Ti)*alphai;
                alphaAlphah[celli] +=
                    t1.kappa(rhoi, ei, Ti)/t1.Cp(rhoi, ei, Ti)*alphai;
                Gamma = t1.Gamma(rhoi, ei, Ti);
                pi = t1.p(rhoi, ei, Ti);
            }
            else if (x1 < this->residualActivation_)
            {
                alphaCp[celli] += t2.Cp(rhoi, ei, Ti)*alphai;
                alphaCv[celli] += t2.Cv(rhoi, ei, Ti)*alphai;
                alphaMu[celli] += t2.mu(rhoi, ei, Ti)*alphai;
                alphaAlphah[celli] +=
                    t2.kappa(rhoi, ei, Ti)/t2.Cp(rhoi, ei, Ti)*alphai;

                Gamma = t2.Gamma(rhoi, ei, Ti);
                pi = t2.p(rhoi, ei, Ti);
            }
            else
            {
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

                Gamma =
                    t1.Gamma(rhoi, ei, Ti)*x1 + t1.Gamma(rhoi, ei, Ti)*x2;
                pi = t1.p(rhoi, ei, Ti)*x1 + t2.p(rhoi, ei, Ti)*x2;
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
        const scalarField px(this->x(patchi));

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
            if (alphai > this->residualAlpha_.value())
            {
                const scalar x2 = px[facei];
                const scalar x1 = 1.0 - x2;
                const scalar& rhoi(prho[facei]);
                const scalar& ei(phe[facei]);
                const scalar& Ti(pT[facei]);
                scalar Gamma;
                scalar pi;

                if (x2 < this->residualActivation_)
                {
                    palphaCp[facei] += t1.Cp(rhoi, ei, Ti)*alphai;
                    palphaCv[facei] += t1.Cv(rhoi, ei, Ti)*alphai;
                    palphaMu[facei] += t1.mu(rhoi, ei, Ti)*alphai;
                    palphaAlphah[facei] +=
                        t1.kappa(rhoi, ei, Ti)/t1.Cp(rhoi, ei, Ti)*alphai;

                    Gamma = t1.Gamma(rhoi, ei, Ti);
                    pi = t1.p(rhoi, ei, Ti);
                }
                else if (x1 < this->residualActivation_)
                {
                    palphaCp[facei] += t2.Cp(rhoi, ei, Ti)*alphai;
                    palphaCv[facei] += t2.Cv(rhoi, ei, Ti)*alphai;
                    palphaMu[facei] += t2.mu(rhoi, ei, Ti)*alphai;
                    palphaAlphah[facei] +=
                        t2.kappa(rhoi, ei, Ti)/t2.Cp(rhoi, ei, Ti)*alphai;

                    Gamma = t2.Gamma(rhoi, ei, Ti);
                    pi = t2.p(rhoi, ei, Ti);
                }
                else
                {
                    palphaCp[facei] +=
                        (
                            t1.Cp(rhoi, ei, Ti)*x1
                          + t2.Cp(rhoi, ei, Ti)*x2
                        )*alphai;
                    palphaCv[facei] +=
                        (
                            t1.Cv(rhoi, ei, Ti)*x1
                          + t2.Cv(rhoi, ei, Ti)*x2
                        )*alphai;
                    palphaMu[facei] +=
                        (
                            t1.mu(rhoi, ei, Ti)*x1
                          + t2.mu(rhoi, ei, Ti)*x2
                        )*alphai;
                    palphaAlphah[facei] +=
                        (
                            t1.kappa(rhoi, ei, Ti)/t1.Cp(rhoi, ei, Ti)*x1
                          + t2.kappa(rhoi, ei, Ti)/t2.Cp(rhoi, ei, Ti)*x2
                        )*alphai;

                    Gamma =
                        t1.Gamma(rhoi, ei, Ti)*x1 + t2.Gamma(rhoi, ei, Ti)*x2;
                    pi = t1.p(rhoi, ei, Ti)*x1 + t2.p(rhoi, ei, Ti)*x2;
                }
                scalar Xii = alphai/(Gamma - 1.0);
                ppXiSum[facei] += pi*Xii;
                pxiSum[facei] += Xii;
            }
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
        const scalar alphai = alpha[celli];
        if (alphai > this->residualAlpha_.value())
        {
            const scalar x2 = this->cellx(celli);
            const scalar x1 = 1.0 - x2;
            const scalar pi = this->p_[celli];
            const scalar rhoi = this->rho_[celli];
            const scalar ei = this->e_[celli];
            const scalar Ti = this->T_[celli];
            scalar cSqr;
            scalar Gamma;

            if (x2 < this->residualActivation_)
            {
                cSqr = t1.cSqr(pi, rhoi, ei, Ti);
                Gamma = t1.Gamma(rhoi, ei, Ti);
            }
            else if (x1 < this->residualActivation_)
            {
                cSqr = t2.cSqr(pi, rhoi, ei, Ti);
                Gamma = t2.Gamma(rhoi, ei, Ti);
            }
            else
            {
                cSqr =
                    t1.cSqr(pi, rhoi, ei, Ti)*x1
                  + t2.cSqr(pi, rhoi, ei, Ti)*x2;
                Gamma =
                    t1.Gamma(rhoi, ei, Ti)*x1 + t2.Gamma(rhoi, ei, Ti)*x2;
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
        const scalarField px(this->x(patchi));
        fvPatchScalarField& pcSqrRhoXiSum =
            cSqrRhoXiSum.boundaryFieldRef()[patchi];

        forAll(pT, facei)
        {
            const scalar alphai = palpha[facei];
            if (alphai > this->residualAlpha_.value())
            {
                const scalar x2 = px[facei];
                const scalar x1 = 1.0 - x2;
                const scalar pi = pp[facei];
                const scalar rhoi = prho[facei];
                const scalar ei = phe[facei];
                const scalar Ti = pT[facei];
                scalar cSqr;
                scalar Gamma;

                if (x2 < this->residualActivation_)
                {
                    cSqr = t1.cSqr(pi, rhoi, ei, Ti);
                    Gamma = t1.Gamma(rhoi, ei, Ti);
                }
                else if (x1 < this->residualActivation_)
                {
                    cSqr = t2.cSqr(pi, rhoi, ei, Ti);
                    Gamma = t2.Gamma(rhoi, ei, Ti);
                }
                else
                {
                    cSqr =
                        t1.cSqr(pi, rhoi, ei, Ti)*x1
                      + t2.cSqr(pi, rhoi, ei, Ti)*x2;
                    Gamma =
                        t1.Gamma(rhoi, ei, Ti)*x1 + t2.Gamma(rhoi, ei, Ti)*x2;
                }
                pcSqrRhoXiSum[facei] += cSqr*rhoi*alphai/(Gamma - 1.0);
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
    this->rho_ == Thermo::blendedVolScalarFieldProperty
    (
        "rho",
        dimDensity,
        &Thermo::thermoType1::rhoPT,
        &Thermo::thermoType2::rhoPT,
        p,
        this->T_
    );
}


template<class Thermo>
Foam::scalar
Foam::detonatingFluidBlastThermo<Thermo>::cellpRhoT(const label celli) const
{
    const scalar& x = this->cellx(celli);
    const scalar rho = this->rho_[celli];
    const scalar e = this->e_[celli];
    const scalar T = this->T_[celli];
    if (x < this->residualActivation_)
    {
        return Thermo::thermoType1::p(rho, e, T);
    }
    else if ((1.0 - x) < this->residualActivation_)
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
    //- Add detonation energy to initially reacted material
    //  restarts are handled in the activation model
    return volScalarField::New
    (
        "eInit",
        Thermo::blendedVolScalarFieldProperty
        (
            "e",
            dimEnergy/dimMass,
            &Thermo::thermoType1::initializeEnergy,
            &Thermo::thermoType2::initializeEnergy,
            p,
            this->rho_,
            this->e_,
            this->T_
        ) + activation_->initESource()
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidBlastThermo<Thermo>::ESource() const
{
    return volScalarField::New
    (
        "ESource",
        (activation_->ESource() + afterburn_->ESource())*this->rho_
    );
}



// ************************************************************************* //
