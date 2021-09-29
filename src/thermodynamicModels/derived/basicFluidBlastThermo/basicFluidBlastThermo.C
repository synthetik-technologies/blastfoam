/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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

#include "basicFluidBlastThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
void Foam::basicFluidBlastThermo<Thermo>::calculate()
{
    const typename Thermo::thermoType& t(*this);
    forAll(this->rho_, celli)
    {
        const scalar& rhoi(this->rho_[celli]);
        scalar& ei(this->he()[celli]);
        scalar& Ti = this->Thermo::baseThermo::T()[celli];

        // Update temperature
        Ti = t.TRhoE(Ti, rhoi, ei);
        if (Ti < this->TLow_)
        {
            ei = t.Es(rhoi, ei, this->TLow_);
            Ti = this->TLow_;
        }

        scalar pi = t.p(rhoi, ei, Ti);
        scalar Cpi = t.Cp(rhoi, ei, Ti);
        this->Thermo::baseThermo::p()[celli] = pi;
        this->CpRef()[celli] = Cpi;
        this->CvRef()[celli] = t.Cv(rhoi, ei, Ti);
        this->muRef()[celli] = t.mu(rhoi, ei, Ti);
        this->alpha()[celli] = t.kappa(rhoi, ei, Ti)/Cpi;
        this->speedOfSound()[celli] =
            sqrt(max(t.cSqr(pi, rhoi, ei, Ti), small));
    }

    this->Thermo::baseThermo::T().correctBoundaryConditions();
    this->he().correctBoundaryConditions();
    this->Thermo::baseThermo::p().correctBoundaryConditions();

    forAll(this->rho_.boundaryField(), patchi)
    {
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pT =
            this->Thermo::baseThermo::T().boundaryField()[patchi];
        const fvPatchScalarField& phe = this->he().boundaryField()[patchi];
        const fvPatchScalarField& pp =
            this->Thermo::baseThermo::p().boundaryField()[patchi];

        fvPatchScalarField& pCp = this->CpRef().boundaryFieldRef()[patchi];
        fvPatchScalarField& pCv = this->CvRef().boundaryFieldRef()[patchi];
        fvPatchScalarField& pmu = this->muRef().boundaryFieldRef()[patchi];
        fvPatchScalarField& palpha =
            this->alpha().boundaryFieldRef()[patchi];
        fvPatchScalarField& pc =
            this->speedOfSound().boundaryFieldRef()[patchi];

        forAll(prho, facei)
        {
            const scalar rhoi(prho[facei]);
            const scalar ei(phe[facei]);
            const scalar Ti(pT[facei]);

            const scalar Cpi = t.Cp(rhoi, ei, Ti);
            pCp[facei] = Cpi;
            pCv[facei] = t.Cv(rhoi, ei, Ti);
            pmu[facei] = t.mu(rhoi, ei, Ti);
            palpha[facei] = t.kappa(rhoi, ei, Ti)/Cpi;
            pc[facei] =
                sqrt(max(t.cSqr(pp[facei], rhoi, ei, Ti), small));
        }
    }
}


template<class Thermo>
void Foam::basicFluidBlastThermo<Thermo>::calculate
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
    const typename Thermo::thermoType& t(*this);
    forAll(alpha, celli)
    {
        const scalar vfi = alpha[celli];
        if (vfi > this->residualAlpha_.value())
        {
            const scalar alphai(alpha[celli]);
            const scalar rhoi(this->rho_[celli]);
            const scalar ei(he[celli]);
            const scalar Ti(T[celli]);
            const scalar Xii = alphai/(t.Gamma(rhoi, ei, Ti) - 1.0);

            alphaCp[celli] += t.Cp(rhoi, ei, Ti)*alphai;
            alphaCv[celli] += t.Cv(rhoi, ei, Ti)*alphai;
            alphaMu[celli] += t.mu(rhoi, ei, Ti)*alphai;
            alphaAlphah[celli] +=
                t.kappa(rhoi, ei, Ti)/t.Cp(rhoi, ei, Ti)*alphai;
            pXiSum[celli] += t.p(rhoi, ei, Ti)*Xii;
            XiSum[celli] += Xii;
            rhoXiSum[celli] += rhoi*Xii;
        }
    }

    forAll(alpha.boundaryField(), patchi)
    {
        const fvPatchScalarField& palpha = alpha.boundaryField()[patchi];
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pT = T.boundaryField()[patchi];
        const fvPatchScalarField& phe = he.boundaryField()[patchi];

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
            const scalar alphai(palpha[facei]);
            if (alphai > this->residualAlpha_.value())
            {
                const scalar rhoi(prho[facei]);
                const scalar ei(phe[facei]);
                const scalar Ti(pT[facei]);
                const scalar Xii = alphai/(t.Gamma(rhoi, ei, Ti) - 1.0);

                const scalar Cpi = t.Cp(rhoi, ei, Ti);

                ppXiSum[facei] = t.p(rhoi, ei, Ti)*Xii;
                palphaCp[facei] = Cpi*alphai;
                palphaCv[facei] = t.Cv(rhoi, ei, Ti)*alphai;
                palphaMu[facei] = t.mu(rhoi, ei, Ti)*alphai;
                palphaAlphah[facei] = t.kappa(rhoi, ei, Ti)/Cpi*alphai;
                pxiSum[facei] += Xii;
                prhoXiSum[facei] += rhoi*Xii;
            }
        }
    }
}


template<class Thermo>
void Foam::basicFluidBlastThermo<Thermo>::calculateSpeedOfSound
(
    const volScalarField& alpha,
    volScalarField& cSqrRhoXiSum
)
{
    const typename Thermo::thermoType& t(*this);
    forAll(this->rho_, celli)
    {
        const scalar vfi = alpha[celli];
        if (vfi > this->residualAlpha_.value())
        {
            cSqrRhoXiSum[celli] +=
                t.cSqr
                (
                    this->p_[celli],
                    this->rho_[celli],
                    this->e_[celli],
                    this->T_[celli]
                )*this->rho_[celli]*vfi
               /(
                   t.Gamma
                   (
                        this->rho_[celli],
                        this->e_[celli],
                        this->T_[celli]
                    )
                 - 1.0
                );
        }
    }

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& palpha = alpha.boundaryField()[patchi];
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        const fvPatchScalarField& phe = this->e_.boundaryField()[patchi];
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pcSqrRhoXiSum =
            cSqrRhoXiSum.boundaryFieldRef()[patchi];

        forAll(pT, facei)
        {
            if (palpha[facei] > this->residualAlpha_.value())
            {
                pcSqrRhoXiSum[facei] +=
                    t.cSqr(pp[facei], prho[facei], phe[facei], pT[facei])
                   *palpha[facei]*prho[facei]
                   /(t.Gamma(prho[facei], phe[facei], pT[facei]) - 1.0);
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::basicFluidBlastThermo<Thermo>::basicFluidBlastThermo
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
        phaseName,
        masterName
    )
{
    //- Initialize the density using the pressure and temperature
    //  This is only done at the first time step (Not on restart)
    if
    (
        max(this->rho_).value() <= 0
     || (
            dict.lookupOrDefault<Switch>("calculateDensity", false)
         && !this->rho_.time().restart()
        )
    )
    {
        updateRho(Thermo::baseThermo::p());
    }
    this->initializeFields();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::basicFluidBlastThermo<Thermo>::~basicFluidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::basicFluidBlastThermo<Thermo>::correct()
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
void Foam::basicFluidBlastThermo<Thermo>::updateRho(const volScalarField& p)
{
    volScalarField rhoNew
    (
        Thermo::volScalarFieldProperty
        (
            "rho",
            dimDensity,
            &Thermo::thermoType::rhoPT,
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
Foam::tmp<Foam::volScalarField>
Foam::basicFluidBlastThermo<Thermo>::ESource() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "ESource",
                this->rho_.mesh().time().timeName(),
                this->rho_.mesh()
            ),
            this->rho_.mesh(),
            dimensionedScalar("0", dimEnergy/dimTime/dimVolume, 0.0)
        )
    );
}


template<class Thermo>
Foam::scalar
Foam::basicFluidBlastThermo<Thermo>::cellGamma(const label celli) const
{
    return Thermo::thermoType::Gamma
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::scalar
Foam::basicFluidBlastThermo<Thermo>::cellpRhoT(const label celli) const
{
    return Thermo::thermoType::p
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::scalar
Foam::basicFluidBlastThermo<Thermo>::celldpdRho(const label celli) const
{
    return Thermo::thermoType::dpdRho
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::scalar
Foam::basicFluidBlastThermo<Thermo>::celldpde(const label celli) const
{
    return Thermo::thermoType::dpde
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::basicFluidBlastThermo<Thermo>::calce(const volScalarField& p) const
{
    return Thermo::volScalarFieldProperty
    (
        "e",
        dimEnergy/dimMass,
        &Thermo::initializeEnergy,
        p,
        this->rho_,
        this->e_,
        this->T_
    );
}

// ************************************************************************* //
