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
    scalarField& eI = this->heRef().primitiveFieldRef();
    scalarField& TI = this->TRef().primitiveFieldRef();
    scalarField& pI = this->pRef().primitiveFieldRef();
    scalarField& CpI = this->CpRef().primitiveFieldRef();
    scalarField& CvI = this->CvRef().primitiveFieldRef();
    scalarField& muI = this->muRef().primitiveFieldRef();
    scalarField& alphaI = this->alphaRef().primitiveFieldRef();
    scalarField& speedOfSoundI = this->speedOfSoundRef().primitiveFieldRef();

    forAll(this->rho_, celli)
    {
        const scalar& rhoi(this->rho_[celli]);
        scalar& ei(eI[celli]);
        scalar& Ti = TI[celli];

        // Update temperature
        Ti = t.TRhoE(Ti, rhoi, ei);
        if (Ti < this->TLow_)
        {
            ei = t.Es(rhoi, ei, this->TLow_);
            Ti = this->TLow_;
        }

        scalar pi = t.p(rhoi, ei, Ti);
        scalar Cpi = t.Cp(rhoi, ei, Ti);
        pI[celli] = pi;
        CpI[celli] = Cpi;
        CvI[celli] = t.Cv(rhoi, ei, Ti);
        muI[celli] = t.mu(rhoi, ei, Ti);
        alphaI[celli] = t.kappa(rhoi, ei, Ti)/Cpi;
        speedOfSoundI[celli] = sqrt(max(t.cSqr(pi, rhoi, ei, Ti), small));
    }

    this->TRef().correctBoundaryConditions();
    this->heRef().correctBoundaryConditions();
    this->pRef().correctBoundaryConditions();

    volScalarField::Boundary& bCp = this->CpRef().boundaryFieldRef();
    volScalarField::Boundary& bCv = this->CvRef().boundaryFieldRef();
    volScalarField::Boundary& bmu = this->muRef().boundaryFieldRef();
    volScalarField::Boundary& balpha = this->alphaRef().boundaryFieldRef();
    volScalarField::Boundary& bspeedOfSound =
        this->speedOfSoundRef().boundaryFieldRef();

    forAll(this->rho_.boundaryField(), patchi)
    {
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->TRef().boundaryField()[patchi];
        const fvPatchScalarField& phe = this->heRef().boundaryField()[patchi];
        const fvPatchScalarField& pp = this->pRef().boundaryField()[patchi];

        fvPatchScalarField& pCp = bCp[patchi];
        fvPatchScalarField& pCv = bCv[patchi];
        fvPatchScalarField& pmu = bmu[patchi];
        fvPatchScalarField& palpha = balpha[patchi];
        fvPatchScalarField& pspeedOfSound = bspeedOfSound[patchi];

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
            pspeedOfSound[facei] =
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
    volScalarField& XiSum
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
        }
    }

    volScalarField::Boundary& balphaCp = alphaCp.boundaryFieldRef();
    volScalarField::Boundary& balphaCv = alphaCv.boundaryFieldRef();
    volScalarField::Boundary& balphaMu = alphaMu.boundaryFieldRef();
    volScalarField::Boundary& balphaAlphah = alphaAlphah.boundaryFieldRef();
    volScalarField::Boundary& bpXiSum = pXiSum.boundaryFieldRef();
    volScalarField::Boundary& bxiSum = XiSum.boundaryFieldRef();

    forAll(alpha.boundaryField(), patchi)
    {
        const fvPatchScalarField& palpha = alpha.boundaryField()[patchi];
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pT = T.boundaryField()[patchi];
        const fvPatchScalarField& phe = he.boundaryField()[patchi];

        fvPatchScalarField& palphaCp = balphaCp[patchi];
        fvPatchScalarField& palphaCv = balphaCv[patchi];
        fvPatchScalarField& palphaMu = balphaMu[patchi];
        fvPatchScalarField& palphaAlphah = balphaAlphah[patchi];
        fvPatchScalarField& ppXiSum = bpXiSum[patchi];
        fvPatchScalarField& pxiSum = bxiSum[patchi];

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
        const scalar alphai = alpha[celli];
        if (alphai > this->residualAlpha_.value())
        {
            cSqrRhoXiSum[celli] +=
                t.cSqr
                (
                    this->p_[celli],
                    this->rho_[celli],
                    this->e_[celli],
                    this->T_[celli]
                )*this->rho_[celli]*alphai
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

    volScalarField::Boundary& bcSqrRhoXiSum = cSqrRhoXiSum.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& palpha = alpha.boundaryField()[patchi];
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        const fvPatchScalarField& phe = this->e_.boundaryField()[patchi];
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        fvPatchScalarField& pcSqrRhoXiSum = bcSqrRhoXiSum[patchi];

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
    this->rho_ == Thermo::volScalarFieldProperty
    (
        "rho",
        dimDensity,
        &Thermo::thermoType::rhoPT,
        this->rho_,
        p,
        this->T_
    );
}


template<class Thermo>
void Foam::basicFluidBlastThermo<Thermo>::updateRho
(
    const volScalarField& alpha,
    const volScalarField& p
)
{
    const typename Thermo::thermoType& t(*this);
    scalarField& rhoI = this->rho_.primitiveFieldRef();
    forAll(rhoI, celli)
    {
        if (alpha[celli] > this->residualAlpha_.value())
        {
            rhoI[celli] = t.rhoPT(rhoI[celli], p[celli], this->T_[celli]);
        }
    }

    volScalarField::Boundary& brho = this->rho_.boundaryFieldRef();

    forAll(this->rho_.boundaryField(), patchi)
    {
        scalarField& prho = brho[patchi];
        const scalarField& palpha = alpha.boundaryField()[patchi];
        const scalarField& pT = this->T_.boundaryField()[patchi];
        const scalarField& pp = p.boundaryField()[patchi];

        forAll(prho, facei)
        {
            if (palpha[facei] > this->residualAlpha_.value())
            {
                prho[facei] = t.rhoPT(prho[facei], pp[facei], pT[facei]);
            }
        }
    }
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::basicFluidBlastThermo<Thermo>::ESource() const
{
    return tmp<volScalarField>
    (
        volScalarField::New
        (
            "ESource",
            this->rho_.mesh(),
            dimensionedScalar("0", dimEnergy/dimTime/dimVolume, 0.0)
        )
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::basicFluidBlastThermo<Thermo>::initESource() const
{
    return tmp<volScalarField>
    (
        volScalarField::New
        (
            "initESource",
            this->rho_.mesh(),
            dimensionedScalar("0", dimEnergy/dimMass, 0.0)
        )
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::basicFluidBlastThermo<Thermo>::Gamma() const
{
    return Thermo::volScalarFieldProperty
    (
        "Gamma",
        dimless,
        &Thermo::thermoType::Gamma,
        this->rho_,
        this->e_,
        this->T_
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
Foam::scalar Foam::basicFluidBlastThermo<Thermo>::cellpRhoT
(
    const label celli,
    const bool limit
) const
{
    return Thermo::thermoType::p
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli],
        limit
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
        &Thermo::thermoType::initializeEnergy,
        p,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class Thermo>
Foam::scalar
Foam::basicFluidBlastThermo<Thermo>::calcCelle
(
    const scalar p,
    const label celli
) const
{
    return Thermo::thermoType::initializeEnergy
    (
        p,
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::basicFluidBlastThermo<Thermo>::calcp() const
{
    return Thermo::volScalarFieldProperty
    (
        "p",
        dimPressure,
        &Thermo::thermoType::pRhoT,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::basicFluidBlastThermo<Thermo>::calcSpeedOfSound() const
{
    tmp<volScalarField> tcSqr
    (
        Thermo::volScalarFieldProperty
        (
            "cSqr",
            sqr(dimVelocity),
            &Thermo::thermoType::cSqr,
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
