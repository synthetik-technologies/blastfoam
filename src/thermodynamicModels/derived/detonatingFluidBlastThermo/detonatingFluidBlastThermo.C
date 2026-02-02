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

    const scalarField& rhoCells = this->rho_.primitiveField();
    scalarField& heCells = this->heRef();
    scalarField& TCells = this->TRef().primitiveFieldRef();
    scalarField& pCells = this->pRef().primitiveFieldRef();
    scalarField& CpCells = this->CpRef().primitiveFieldRef();
    scalarField& CvCells = this->CvRef().primitiveFieldRef();
    scalarField& muCells = this->muRef().primitiveFieldRef();
    scalarField& kappaCells = this->kappaRef().primitiveFieldRef();
    scalarField& speedOfSoundCells =
        this->speedOfSoundRef().primitiveFieldRef();

    forAll(this->rho_, celli)
    {
        const scalar x2 = this->cellx(celli);
        const scalar x1 = 1.0 - x2;
        const scalar rhoi = rhoCells[celli];
        scalar& ei = heCells[celli];
        scalar& Ti = TCells[celli];

        if (x2 < this->residualFac_)
        {
            const scalar eTest = t1.Es(rhoi, ei, this->TLow_);
            if (ei < eTest)
            {
                ei = eTest;
                Ti = this->TLow_;
            }
            else
            {
                Ti = t1.TRhoE(Ti, rhoi, ei);
            }

            const scalar pi = t1.p(rhoi, ei, Ti);

            pCells[celli] = pi;
            CpCells[celli] = t1.Cp(rhoi, ei, Ti);
            CvCells[celli] = t1.Cv(rhoi, ei, Ti);
            muCells[celli] = t1.mu(rhoi, ei, Ti);
            kappaCells[celli] = t1.kappa(rhoi, ei, Ti);
            speedOfSoundCells[celli] =
                sqrt(max(t1.cSqr(pi, rhoi, ei, Ti), small));
        }
        else if (x1 < this->residualFac_)
        {
            const scalar eTest = t2.Es(rhoi, ei, this->TLow_);
            if (ei < eTest)
            {
                ei = eTest;
                Ti = this->TLow_;
            }
            else
            {
                Ti = t2.TRhoE(Ti, rhoi, ei);
            }

            const scalar pi = t2.p(rhoi, ei, Ti);

            pCells[celli] = pi;
            CpCells[celli] = t2.Cp(rhoi, ei, Ti);
            CvCells[celli] = t2.Cv(rhoi, ei, Ti);
            muCells[celli] = t2.mu(rhoi, ei, Ti);
            kappaCells[celli] = t2.kappa(rhoi, ei, Ti);
            speedOfSoundCells[celli] =
                sqrt(max(t2.cSqr(pi, rhoi, ei, Ti), small));
        }
        else
        {

            const scalar eTest =
                t1.Es(rhoi, ei, this->TLow_)*x1
              + t2.Es(rhoi, ei, this->TLow_)*x2;
            if (ei < eTest)
            {
                ei = eTest;
                Ti = this->TLow_;
            }
            else
            {
                Ti =
                    t1.TRhoE(Ti, rhoi, ei)*x1
                  + t2.TRhoE(Ti, rhoi, ei)*x2;
            }

            const scalar pi =
                t1.p(rhoi, ei, Ti)*x1
              + t2.p(rhoi, ei, Ti)*x2;

            pCells[celli] = pi;
            CpCells[celli] =
                t1.Cp(rhoi, ei, Ti)*x1
              + t2.Cp(rhoi, ei, Ti)*x2;
            CvCells[celli] =
                t1.Cv(rhoi, ei, Ti)*x1
              + t2.Cv(rhoi, ei, Ti)*x2;
            muCells[celli] =
                t1.mu(rhoi, ei, Ti)*x1
              + t2.mu(rhoi, ei, Ti)*x2;
            kappaCells[celli] =
                t1.kappa(rhoi, ei, Ti)*x1
              + t2.kappa(rhoi, ei, Ti)*x2;
            speedOfSoundCells[celli] =
                sqrt
                (
                    max(t1.cSqr(pi, rhoi, ei, Ti), small)*x1
                  + max(t2.cSqr(pi, rhoi, ei, Ti), small)*x2
                );
        }
    }

    const volScalarField::Boundary& rhoBf = this->rho_.boundaryField();

    volScalarField::Boundary& heBf = this->heRef().boundaryFieldRef();
    volScalarField::Boundary& TBf = this->TRef().boundaryFieldRef();
    volScalarField::Boundary& pBf = this->pRef().boundaryFieldRef();

    volScalarField::Boundary& CpBf = this->CpRef().boundaryFieldRef();
    volScalarField::Boundary& CvBf = this->CvRef().boundaryFieldRef();
    volScalarField::Boundary& muBf = this->muRef().boundaryFieldRef();
    volScalarField::Boundary& kappaBf =
        this->kappaRef().boundaryFieldRef();
    volScalarField::Boundary& speedOfSoundBf =
        this->speedOfSoundRef().boundaryFieldRef();

    this->pRef().correctBoundaryConditions();

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& prho = rhoBf[patchi];
        const fvPatchScalarField& pp = pBf[patchi];
        tmp<scalarField> tpx(this->x(patchi));
        const scalarField& px = tpx();

        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];
        fvPatchScalarField& pc = speedOfSoundBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const scalar x2 = px[facei];
                const scalar x1 = 1.0 - x2;
                const scalar rhoi = prho[facei];
                const scalar pi = pp[facei];

                scalar& ei = phe[facei];
                scalar& Ti = pT[facei];

                if (x2 < this->residualFac_)
                {
                    phe[facei] = t1.Es(rhoi, ei, Ti);
                    pCp[facei] = t1.Cp(rhoi, ei, Ti);
                    pCv[facei] = t1.Cv(rhoi, ei, Ti);
                    pmu[facei] = t1.mu(rhoi, ei, Ti);
                    pkappa[facei] = t1.kappa(rhoi, ei, Ti);
                    pc[facei] = sqrt(max(t1.cSqr(pi, rhoi, ei, Ti), small));
                }
                else if (x1 < this->residualFac_)
                {
                    phe[facei] = t2.Es(rhoi, ei, Ti);
                    pCp[facei] = t2.Cp(rhoi, ei, Ti);
                    pCv[facei] = t2.Cv(rhoi, ei, Ti);
                    pmu[facei] = t2.mu(rhoi, ei, Ti);
                    pkappa[facei] = t2.kappa(rhoi, ei, Ti);
                    pc[facei] = sqrt(max(t2.cSqr(pi, rhoi, ei, Ti), small));
                }
                else
                {
                    phe[facei] =
                        t1.Es(rhoi, ei, Ti)*x1
                      + t2.Es(rhoi, ei, Ti)*x2;
                    pCp[facei] =
                        t1.Cp(rhoi, ei, Ti)*x1
                      + t2.Cp(rhoi, ei, Ti)*x2;
                    pCv[facei] =
                        t1.Cv(rhoi, ei, Ti)*x1
                      + t2.Cv(rhoi, ei, Ti)*x2;
                    pmu[facei] =
                        t1.mu(rhoi, ei, Ti)*x1
                      + t2.mu(rhoi, ei, Ti)*x2;
                    pkappa[facei] =
                        t1.kappa(rhoi, ei, Ti)*x1
                      + t2.kappa(rhoi, ei, Ti)*x2;
                    pc[facei] =
                        sqrt
                        (
                            max(t1.cSqr(pi, rhoi, ei, Ti), small)*x1
                          + max(t2.cSqr(pi, rhoi, ei, Ti), small)*x2
                        );
                }
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const scalar x2 = px[facei];
                const scalar x1 = 1.0 - x2;

                const scalar rhoi = prho[facei];
                const scalar pi = pp[facei];

                scalar& ei = phe[facei];
                scalar& Ti = pT[facei];

                if (x2 < this->residualFac_)
                {
                    const scalar eTest = t1.Es(rhoi, ei, Ti);
                    if (Ti < this->TLow_)
                    {
                        Ti = this->TLow_;
                        ei = eTest;
                    }
                    else
                    {
                        Ti = t1.TRhoE(Ti, rhoi, ei);
                    }
                    pCp[facei] = t1.Cp(rhoi, ei, Ti);
                    pCv[facei] = t1.Cv(rhoi, ei, Ti);
                    pmu[facei] = t1.mu(rhoi, ei, Ti);
                    pkappa[facei] = t1.kappa(rhoi, ei, Ti);
                    pc[facei] = sqrt(max(t1.cSqr(pi, rhoi, ei, Ti), small));
                }
                else if (x1 < this->residualFac_)
                {
                    const scalar eTest = t2.Es(rhoi, ei, Ti);
                    if (Ti < this->TLow_)
                    {
                        Ti = this->TLow_;
                        ei = eTest;
                    }
                    else
                    {
                        Ti = t2.TRhoE(Ti, rhoi, ei);
                    }
                    pCp[facei] = t2.Cp(rhoi, ei, Ti);
                    pCv[facei] = t2.Cv(rhoi, ei, Ti);
                    pmu[facei] = t2.mu(rhoi, ei, Ti);
                    pkappa[facei] = t2.kappa(rhoi, ei, Ti);
                    pc[facei] = sqrt(max(t2.cSqr(pi, rhoi, ei, Ti), small));
                }
                else
                {
                    const scalar eTest =
                        t1.Es(rhoi, ei, this->TLow_)*x1
                      + t2.Es(rhoi, ei, this->TLow_)*x2;
                    if (ei < eTest)
                    {
                        ei = eTest;
                        Ti = this->TLow_;
                    }
                    else
                    {
                        Ti =
                            t1.TRhoE(Ti, rhoi, ei)*x1
                          + t2.TRhoE(Ti, rhoi, ei)*x2;
                    }
                    pCp[facei] =
                        t1.Cp(rhoi, ei, Ti)*x1
                      + t2.Cp(rhoi, ei, Ti)*x2;
                    pCv[facei] =
                        t1.Cv(rhoi, ei, Ti)*x1
                      + t2.Cv(rhoi, ei, Ti)*x2;
                    pmu[facei] =
                        t1.mu(rhoi, ei, Ti)*x1
                      + t2.mu(rhoi, ei, Ti)*x2;
                    pkappa[facei] =
                        t1.kappa(rhoi, ei, Ti)*x1
                      + t2.kappa(rhoi, ei, Ti)*x2;
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
    volScalarField& alphaKappa,
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
        const scalar rhoi = this->rho_[celli];
        const scalar ei = he[celli];
        const scalar Ti = T[celli];
        if (alphai > this->residualAlpha_.value())
        {
            scalar Gamma, pi;
            if (x2 < this->residualFac_)
            {
                alphaCp[celli] += t1.Cp(rhoi, ei, Ti)*alphai*rhoi;
                alphaCv[celli] += t1.Cv(rhoi, ei, Ti)*alphai*rhoi;
                alphaMu[celli] += t1.mu(rhoi, ei, Ti)*alphai*rhoi;
                alphaKappa[celli] += t1.kappa(rhoi, ei, Ti)*alphai*rhoi;
                Gamma = t1.Gamma(rhoi, ei, Ti);
                pi = t1.p(rhoi, ei, Ti);
            }
            else if (x1 < this->residualFac_)
            {
                alphaCp[celli] += t2.Cp(rhoi, ei, Ti)*alphai*rhoi;
                alphaCv[celli] += t2.Cv(rhoi, ei, Ti)*alphai*rhoi;
                alphaMu[celli] += t2.mu(rhoi, ei, Ti)*alphai*rhoi;
                alphaKappa[celli] += t2.kappa(rhoi, ei, Ti)*alphai*rhoi;

                Gamma = t2.Gamma(rhoi, ei, Ti);
                pi = t2.p(rhoi, ei, Ti);
            }
            else
            {
                alphaCp[celli] +=
                    (
                        t1.Cp(rhoi, ei, Ti)*x1
                      + t2.Cp(rhoi, ei, Ti)*x2
                    )*alphai*rhoi;
                alphaCv[celli] +=
                    (
                        t1.Cv(rhoi, ei, Ti)*x1
                      + t2.Cv(rhoi, ei, Ti)*x2
                    )*alphai*rhoi;
                alphaMu[celli] +=
                    (
                        t1.mu(rhoi, ei, Ti)*x1
                      + t2.mu(rhoi, ei, Ti)*x2
                    )*alphai*rhoi;
                alphaKappa[celli] +=
                    (
                        t1.kappa(rhoi, ei, Ti)*x1
                      + t2.kappa(rhoi, ei, Ti)*x2
                    )*alphai*rhoi;

                Gamma =
                    t1.Gamma(rhoi, ei, Ti)*x1 + t1.Gamma(rhoi, ei, Ti)*x2;
                pi = t1.p(rhoi, ei, Ti)*x1 + t2.p(rhoi, ei, Ti)*x2;
            }
            scalar Xii = alphai/Gamma;
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
        fvPatchScalarField& palphaKappa =
            alphaKappa.boundaryFieldRef()[patchi];
        fvPatchScalarField& ppXiSum = pXiSum.boundaryFieldRef()[patchi];
        fvPatchScalarField& pxiSum = XiSum.boundaryFieldRef()[patchi];

        forAll(palpha, facei)
        {
            const scalar alphai = palpha[facei];
            if (alphai > this->residualAlpha_.value())
            {
                const scalar x2 = px[facei];
                const scalar x1 = 1.0 - x2;

                const scalar rhoi = prho[facei];
                const scalar ei = phe[facei];
                const scalar Ti = pT[facei];

                scalar Gamma, pi;
                if (x2 < this->residualFac_)
                {
                    palphaCp[facei] += t1.Cp(rhoi, ei, Ti)*alphai*rhoi;
                    palphaCv[facei] += t1.Cv(rhoi, ei, Ti)*alphai*rhoi;
                    palphaMu[facei] += t1.mu(rhoi, ei, Ti)*alphai*rhoi;
                    palphaKappa[facei] += t1.kappa(rhoi, ei, Ti)*alphai*rhoi;

                    Gamma = t1.Gamma(rhoi, ei, Ti);
                    pi = t1.p(rhoi, ei, Ti);
                }
                else if (x1 < this->residualFac_)
                {
                    palphaCp[facei] += t2.Cp(rhoi, ei, Ti)*alphai*rhoi;
                    palphaCv[facei] += t2.Cv(rhoi, ei, Ti)*alphai*rhoi;
                    palphaMu[facei] += t2.mu(rhoi, ei, Ti)*alphai*rhoi;
                    palphaKappa[facei] += t2.kappa(rhoi, ei, Ti)*alphai*rhoi;

                    Gamma = t2.Gamma(rhoi, ei, Ti);
                    pi = t2.p(rhoi, ei, Ti);
                }
                else
                {
                    palphaCp[facei] +=
                        (
                            t1.Cp(rhoi, ei, Ti)*x1
                          + t2.Cp(rhoi, ei, Ti)*x2
                        )*alphai*rhoi;
                    palphaCv[facei] +=
                        (
                            t1.Cv(rhoi, ei, Ti)*x1
                          + t2.Cv(rhoi, ei, Ti)*x2
                        )*alphai*rhoi;
                    palphaMu[facei] +=
                        (
                            t1.mu(rhoi, ei, Ti)*x1
                          + t2.mu(rhoi, ei, Ti)*x2
                        )*alphai*rhoi;
                    palphaKappa[facei] +=
                        (
                            t1.kappa(rhoi, ei, Ti)*x1
                          + t2.kappa(rhoi, ei, Ti)*x2
                        )*alphai*rhoi;

                    Gamma =
                        t1.Gamma(rhoi, ei, Ti)*x1 + t2.Gamma(rhoi, ei, Ti)*x2;
                    pi = t1.p(rhoi, ei, Ti)*x1 + t2.p(rhoi, ei, Ti)*x2;
                }
                scalar Xii = alphai/Gamma;
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

            scalar cSqr, Gamma;
            if (x2 < this->residualFac_)
            {
                cSqr = t1.cSqr(pi, rhoi, ei, Ti);
                Gamma = t1.Gamma(rhoi, ei, Ti);
            }
            else if (x1 < this->residualFac_)
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
            cSqrRhoXiSum[celli] += cSqr*rhoi*alphai/Gamma;
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

                scalar cSqr, Gamma;
                if (x2 < this->residualFac_)
                {
                    cSqr = t1.cSqr(pi, rhoi, ei, Ti);
                    Gamma = t1.Gamma(rhoi, ei, Ti);
                }
                else if (x1 < this->residualFac_)
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
                pcSqrRhoXiSum[facei] += cSqr*rhoi*alphai/Gamma;
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
    dict.readIfPresent("residualActivation", this->residualFac_);

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
void Foam::detonatingFluidBlastThermo<Thermo>::solveExplicit()
{
    activation_->solveExplicit();
    afterburn_->solveExplicit();
}


template<class Thermo>
void Foam::detonatingFluidBlastThermo<Thermo>::storeExplicit()
{
    activation_->storeExplicit();
    afterburn_->storeExplicit();
}


template<class Thermo>
void Foam::detonatingFluidBlastThermo<Thermo>::solveImplicit()
{
    activation_->solveImplicit();
    afterburn_->solveImplicit();
}


template<class Thermo>
void Foam::detonatingFluidBlastThermo<Thermo>::clear()
{
    activation_->clear();
    afterburn_->clear();
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
        this->rho_,
        p,
        this->T_
    );
}


template<class Thermo>
void Foam::detonatingFluidBlastThermo<Thermo>::updateRho
(
    const volScalarField& alpha,
    const volScalarField& p
)
{
    const typename Thermo::thermoType1& t1(*this);
    const typename Thermo::thermoType2& t2(*this);

    scalarField& rhoI = this->rho_.primitiveFieldRef();
    forAll(this->rho_, celli)
    {
        if (alpha[celli] > this->residualAlpha_.value())
        {
            const scalar x2 = this->cellx(celli);
            const scalar x1 = 1.0 - x2;

            if (x2 < this->residualFac_)
            {
                rhoI[celli] = t1.rhoPT(rhoI[celli], p[celli], this->T_[celli]);
            }
            else if (x1 < this->residualFac_)
            {
                rhoI[celli] = t2.rhoPT(rhoI[celli], p[celli], this->T_[celli]);
            }
            else
            {
                rhoI[celli] =
                    t1.rhoPT(rhoI[celli], p[celli], this->T_[celli])*x1
                  + t2.rhoPT(rhoI[celli], p[celli], this->T_[celli])*x2;
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
        const scalarField px(this->x(patchi));

        forAll(prho, facei)
        {
            if (palpha[facei] > this->residualAlpha_.value())
            {
                const scalar x2 = px[facei];
                const scalar x1 = 1.0 - x2;
                if (x2 < this->residualFac_)
                {
                    prho[facei] = t1.rhoPT(prho[facei], pp[facei], pT[facei]);
                }
                else if (x1 < this->residualFac_)
                {
                    prho[facei] = t2.rhoPT(prho[facei], pp[facei], pT[facei]);
                }
                else
                {
                    prho[facei] =
                        t1.rhoPT(prho[facei], pp[facei], pT[facei])*x1
                      + t2.rhoPT(prho[facei], pp[facei], pT[facei])*x2;
                }
            }
        }
    }
}


template<class Thermo>
Foam::scalar Foam::detonatingFluidBlastThermo<Thermo>::cellpRhoT
(
    const label celli,
    const bool limit
) const
{
    const scalar& x = this->cellx(celli);
    const scalar rho = this->rho_[celli];
    const scalar e = this->e_[celli];
    const scalar T = this->T_[celli];
    if (x < this->residualFac_)
    {
        return Thermo::thermoType1::p(rho, e, T, limit);
    }
    else if ((1.0 - x) < this->residualFac_)
    {
        return Thermo::thermoType2::p(rho, e, T, limit);
    }

    return
        Thermo::thermoType2::p(rho, e, T, limit)*x
      + Thermo::thermoType1::p(rho, e, T, limit)*(1.0 - x);
}


template<class Thermo>
Foam::scalar Foam::detonatingFluidBlastThermo<Thermo>::patchFacepRhoT
(
    const label patchi,
    const label facei,
    const bool limit
) const
{
    const scalar& x = this->patchFacex(patchi, facei);
    const scalar rho = this->rho_.boundaryField()[patchi][facei];
    const scalar e = this->e_.boundaryField()[patchi][facei];
    const scalar T = this->T_.boundaryField()[patchi][facei];
    if (x < this->residualFac_)
    {
        return Thermo::thermoType1::p(rho, e, T, limit);
    }
    else if ((1.0 - x) < this->residualFac_)
    {
        return Thermo::thermoType2::p(rho, e, T, limit);
    }

    return
        Thermo::thermoType2::p(rho, e, T, limit)*x
      + Thermo::thermoType1::p(rho, e, T, limit)*(1.0 - x);
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidBlastThermo<Thermo>::Gamma() const
{
    return Thermo::blendedVolScalarFieldProperty
    (
        "Gamma",
        dimless,
        &Thermo::thermoType1::Gamma,
        &Thermo::thermoType2::Gamma,
        this->rho_,
        this->e_,
        this->T_
    );
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
Foam::detonatingFluidBlastThermo<Thermo>::patchFaceGamma
(
    const label patchi,
    const label facei
) const
{
    const scalar& x = this->patchFacex(patchi, facei);
    const scalar rho = this->rho_.boundaryField()[patchi][facei];
    const scalar e = this->e_.boundaryField()[patchi][facei];
    const scalar T = this->T_.boundaryField()[patchi][facei];
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
Foam::detonatingFluidBlastThermo<Thermo>::cellSpeedOfSound
(
    const scalar p,
    const label celli
) const
{
    const scalar& x = this->cellx(celli);
    const scalar rho = this->rho_[celli];
    const scalar e = this->e_[celli];
    const scalar T = this->T_[celli];
    if (x < small)
    {
        return sqrt(Thermo::thermoType1::cSqr(p, rho, e, T));
    }
    else if ((1.0 - x) < small)
    {
        return sqrt(Thermo::thermoType2::cSqr(p, rho, e, T));
    }

    return sqrt
    (
        Thermo::thermoType2::cSqr(p, rho, e, T)*x
      + Thermo::thermoType1::cSqr(p, rho, e, T)*(1.0 - x)
    );
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
Foam::scalar
Foam::detonatingFluidBlastThermo<Thermo>::celldpdT(const label celli) const
{
    const scalar& x = this->cellx(celli);
    const scalar rho = this->rho_[celli];
    const scalar e = this->e_[celli];
    const scalar T = this->T_[celli];
    if (x < small)
    {
        return Thermo::thermoType1::dpdT(rho, e, T);
    }
    else if ((1.0 - x) < small)
    {
        return Thermo::thermoType2::dpdT(rho, e, T);
    }

    return
        Thermo::thermoType2::dpdT(rho, e, T)*x
      + Thermo::thermoType1::dpdT(rho, e, T)*(1.0 - x);
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
Foam::scalar
Foam::detonatingFluidBlastThermo<Thermo>::calcCelle
(
    const scalar p,
    const label celli
) const
{
    const scalar& x = this->cellx(celli);
    const scalar rho = this->rho_[celli];
    const scalar e = this->e_[celli];
    const scalar T = this->T_[celli];
    if (x < small)
    {
        return Thermo::thermoType1::initializeEnergy(p, rho, e, T);
    }
    else if ((1.0 - x) < small)
    {
        return Thermo::thermoType2::initializeEnergy(p, rho, e, T);
    }

    return
        Thermo::thermoType2::initializeEnergy(p, rho, e, T)*x
      + Thermo::thermoType1::initializeEnergy(p, rho, e, T)*(1.0 - x);
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


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidBlastThermo<Thermo>::ESource
(
    const volScalarField& alpha
) const
{
    return alpha*ESource();
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidBlastThermo<Thermo>::initESource() const
{
    return activation_->initESource();
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidBlastThermo<Thermo>::initESource
(
    const volScalarField& alpha
) const
{
    return alpha*initESource();
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidBlastThermo<Thermo>::calcp() const
{
    return Thermo::blendedVolScalarFieldProperty
    (
        "p",
        dimPressure,
        &Thermo::thermoType1::pRhoT,
        &Thermo::thermoType2::pRhoT,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidBlastThermo<Thermo>::calcSpeedOfSound() const
{
    tmp<volScalarField> tcSqr
    (
        Thermo::blendedVolScalarFieldProperty
        (
            "cSqr",
            sqr(dimVelocity),
            &Thermo::thermoType1::cSqr,
            &Thermo::thermoType2::cSqr,
            this->p_,
            this->rho_,
            this->e_,
            this->T_
        )
    );
    tcSqr.ref().max(small);
    return sqrt(tcSqr);
}


// ************************************************************************* //
