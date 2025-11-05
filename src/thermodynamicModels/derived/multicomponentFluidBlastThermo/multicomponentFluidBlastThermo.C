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

#include "multicomponentFluidBlastThermo.H"
#include "fvc.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
void Foam::multicomponentFluidBlastThermo<Thermo>::calculate()
{
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
        const typename Thermo::thermoType& t =
            this->cellMixture(celli);

        const scalar& rhoi = rhoCells[celli];
        scalar& ei = heCells[celli];
        scalar& Ti = TCells[celli];

        // Update temperature
        Ti = t.TRhoE(Ti, rhoi, ei);
        if (Ti < this->TLow_)
        {
            ei = t.Es(rhoi, ei, this->TLow_);
            Ti = this->TLow_;
        }

        const scalar pi = t.p(rhoi, ei, Ti);
        pCells[celli] = pi;
        CpCells[celli] = t.Cp(rhoi, ei, Ti);
        CvCells[celli] = t.Cv(rhoi, ei, Ti);
        muCells[celli] = t.mu(rhoi, ei, Ti);
        kappaCells[celli] = t.kappa(rhoi, ei, Ti);
        speedOfSoundCells[celli] =
            sqrt(max(t.cSqr(pi, rhoi, ei, Ti), small));
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

    forAll(this->rho_.boundaryField(), patchi)
    {
        const fvPatchScalarField& prho = rhoBf[patchi];
        const fvPatchScalarField& pp = pBf[patchi];

        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];
        fvPatchScalarField& pc = speedOfSoundBf[patchi];

        if (pT.fixesValue())
        {
            forAll(prho, facei)
            {
                const typename Thermo::thermoType& t =
                    this->patchFaceMixture(patchi, facei);

                const scalar rhoi = prho[facei];
                const scalar Ti = pT[facei];
                scalar& ei = phe[facei];

                ei = t.Es(rhoi, ei, Ti);
                pCp[facei] = t.Cp(rhoi, ei, Ti);
                pCv[facei] = t.Cv(rhoi, ei, Ti);
                pmu[facei] = t.mu(rhoi, ei, Ti);
                pkappa[facei] = t.kappa(rhoi, ei, Ti);
                pc[facei] =
                    sqrt(max(t.cSqr(pp[facei], rhoi, ei, Ti), small));
            }
        }
        else
        {
            forAll(prho, facei)
            {
                const typename Thermo::thermoType& t =
                    this->patchFaceMixture(patchi, facei);

                const scalar rhoi = prho[facei];
                scalar& ei = phe[facei];
                scalar& Ti = pT[facei];

                Ti = t.TRhoE(Ti, rhoi, ei);
                if (Ti < this->TLow_)
                {
                    Ti = this->TLow_;
                    ei = t.Es(rhoi, ei, Ti);
                }
                pCp[facei] = t.Cp(rhoi, ei, Ti);
                pCv[facei] = t.Cv(rhoi, ei, Ti);
                pmu[facei] = t.mu(rhoi, ei, Ti);
                pkappa[facei] = t.kappa(rhoi, ei, Ti);
                pc[facei] =
                    sqrt(max(t.cSqr(pp[facei], rhoi, ei, Ti), small));
            }
        }
    }
}


template<class Thermo>
void Foam::multicomponentFluidBlastThermo<Thermo>::calculate
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
    forAll(alpha, celli)
    {
        const scalar vfi = alpha[celli];
        if (vfi > this->residualAlpha_.value())
        {
            const typename Thermo::thermoType& t =
                this->cellMixture(celli);

            const scalar alphai = alpha[celli];
            const scalar rhoi = this->rho_[celli];
            const scalar ei = he[celli];
            const scalar Ti = T[celli];
            const scalar Xii = alphai/t.Gamma(rhoi, ei, Ti);

            alphaCp[celli] += t.Cp(rhoi, ei, Ti)*alphai;
            alphaCv[celli] += t.Cv(rhoi, ei, Ti)*alphai;
            alphaMu[celli] += t.mu(rhoi, ei, Ti)*alphai;
            alphaKappa[celli] += t.kappa(rhoi, ei, Ti)*alphai;
            pXiSum[celli] += t.p(rhoi, ei, Ti)*Xii;
            XiSum[celli] += Xii;
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
        fvPatchScalarField& palphaKappa =
            alphaKappa.boundaryFieldRef()[patchi];
        fvPatchScalarField& ppXiSum = pXiSum.boundaryFieldRef()[patchi];
        fvPatchScalarField& pxiSum = XiSum.boundaryFieldRef()[patchi];

        forAll(palpha, facei)
        {
            const scalar alphai(palpha[facei]);
            if (alphai > this->residualAlpha_.value())
            {
                const typename Thermo::thermoType& t =
                    this->patchFaceMixture(patchi, facei);

                const scalar rhoi = prho[facei];
                const scalar ei = phe[facei];
                const scalar Ti = pT[facei];
                const scalar Xii = alphai/t.Gamma(rhoi, ei, Ti);

                ppXiSum[facei] = t.p(rhoi, ei, Ti)*Xii;
                palphaCp[facei] = t.Cp(rhoi, ei, Ti)*alphai;
                palphaCv[facei] = t.Cv(rhoi, ei, Ti)*alphai;
                palphaMu[facei] = t.mu(rhoi, ei, Ti)*alphai;
                palphaKappa[facei] = t.kappa(rhoi, ei, Ti)*alphai;
                pxiSum[facei] += Xii;
            }
        }
    }
}


template<class Thermo>
void Foam::multicomponentFluidBlastThermo<Thermo>::calculateSpeedOfSound
(
    const volScalarField& alpha,
    volScalarField& cSqrRhoXiSum
)
{
    forAll(this->rho_, celli)
    {
        const scalar vfi = alpha[celli];
        if (vfi > this->residualAlpha_.value())
        {
            const typename Thermo::thermoType& t =
                this->cellMixture(celli);
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
                const typename Thermo::thermoType& t =
                    this->patchFaceMixture(patchi, facei);

                pcSqrRhoXiSum[facei] +=
                    t.cSqr(pp[facei], prho[facei], phe[facei], pT[facei])
                   *palpha[facei]*prho[facei]
                   /t.Gamma(prho[facei], phe[facei], pT[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::multicomponentFluidBlastThermo<Thermo>::multicomponentFluidBlastThermo
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


template<class Thermo>
Foam::multicomponentFluidBlastThermo<Thermo>::multicomponentFluidBlastThermo
(
    const HashPtrTable<typename Thermo::thermoType, word, string::hash>& thermoData,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
:
    Thermo
    (
        thermoData,
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
Foam::multicomponentFluidBlastThermo<Thermo>::~multicomponentFluidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::multicomponentFluidBlastThermo<Thermo>::correct()
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
void Foam::multicomponentFluidBlastThermo<Thermo>::updateRho
(
    const volScalarField& p
)
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
void Foam::multicomponentFluidBlastThermo<Thermo>::updateRho
(
    const volScalarField& alpha,
    const volScalarField& p
)
{
    scalarField& rhoI = this->rho_.primitiveFieldRef();
    forAll(rhoI, celli)
    {
        if (alpha[celli] > this->residualAlpha_.value())
        {
            const typename Thermo::thermoType& t(this->cellMixture(celli));
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
                const typename Thermo::thermoType& t
                (
                    this->patchFaceMixture(patchi, facei)
                );
                prho[facei] = t.rhoPT(prho[facei], pp[facei], pT[facei]);
            }
        }
    }
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidBlastThermo<Thermo>::ESource() const
{
    return volScalarField::New
    (
        "ESource",
        this->rho_.mesh(),
        dimensionedScalar("0", dimEnergy/dimTime/dimVolume, 0.0)
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidBlastThermo<Thermo>::ESource
(
    const volScalarField& alpha
) const
{
    return ESource();
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidBlastThermo<Thermo>::initESource() const
{
    return volScalarField::New
    (
        "initESource",
        this->rho_.mesh(),
        dimensionedScalar("0", dimEnergy/dimMass, 0.0)
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidBlastThermo<Thermo>::initESource
(
    const volScalarField& alpha
) const
{
    return initESource();
}


template<class Thermo>
bool Foam::multicomponentFluidBlastThermo<Thermo>::inviscid() const
{
    forAll(this->speciesData_, i)
    {
        if (!this->speciesData_[i].inviscid())
        {
            return false;
        }
    }
    return true;
}


template<class Thermo>
Foam::scalar Foam::multicomponentFluidBlastThermo<Thermo>::cellSpeedOfSound
(
    const scalar p,
    const label celli
) const
{
    return sqrt
    (
        this->cellMixture(celli).cSqr
        (
            p,
            this->rho_[celli],
            this->e_[celli],
            this->T_[celli]
        )
    );
}


template<class Thermo>
Foam::scalar Foam::multicomponentFluidBlastThermo<Thermo>::cellpRhoT
(
    const label celli,
    const bool limit
) const
{
    return this->cellMixture(celli).p
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli],
        limit
    );
}


template<class Thermo>
Foam::scalar Foam::multicomponentFluidBlastThermo<Thermo>::patchFacepRhoT
(
    const label patchi,
    const label facei,
    const bool limit
) const
{
    return this->patchFaceMixture(patchi, facei).p
    (
        this->rho_.boundaryField()[patchi][facei],
        this->e_.boundaryField()[patchi][facei],
        this->T_.boundaryField()[patchi][facei],
        limit
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidBlastThermo<Thermo>::Gamma() const
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
Foam::multicomponentFluidBlastThermo<Thermo>::cellGamma(const label celli) const
{
    return this->cellMixture(celli).Gamma
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::scalar Foam::multicomponentFluidBlastThermo<Thermo>::patchFaceGamma
(
    const label patchi,
    const label facei
) const
{
    return this->patchFaceMixture(patchi, facei).Gamma
    (
        this->rho_.boundaryField()[patchi][facei],
        this->e_.boundaryField()[patchi][facei],
        this->T_.boundaryField()[patchi][facei]
    );
}


template<class Thermo>
Foam::scalar
Foam::multicomponentFluidBlastThermo<Thermo>::celldpdRho(const label celli) const
{
    return this->cellMixture(celli).dpdRho
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::scalar
Foam::multicomponentFluidBlastThermo<Thermo>::celldpde(const label celli) const
{
    return this->cellMixture(celli).dpde
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::scalar
Foam::multicomponentFluidBlastThermo<Thermo>::celldpdT(const label celli) const
{
    return this->cellMixture(celli).dpdT
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidBlastThermo<Thermo>::calce
(
    const volScalarField& p
) const
{
    return this->volScalarFieldProperty
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
Foam::multicomponentFluidBlastThermo<Thermo>::calcCelle
(
    const scalar p,
    const label celli
) const
{
    return this->cellMixture(celli).initializeEnergy
    (
        p,
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidBlastThermo<Thermo>::calcp() const
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
Foam::multicomponentFluidBlastThermo<Thermo>::calcSpeedOfSound() const
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


template<class Thermo>
Foam::scalar
Foam::multicomponentFluidBlastThermo<Thermo>::pi
(
    const label speciei,
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->speciesData_[speciei].pRhoT(rho, e, T);
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidBlastThermo<Thermo>::pi
(
    const label speciei,
    const volScalarField& rho,
    const volScalarField& e,
    const volScalarField& T
) const
{
    return this->volScalarFieldSpecieProperty
    (
        speciei,
        "p",
        dimPressure,
        &Thermo::thermoType::pRhoT,
        rho,
        e,
        T
    );
}


template<class Thermo>
Foam::scalar
Foam::multicomponentFluidBlastThermo<Thermo>::dpdRhoi
(
    const label speciei,
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return
      - this->speciesData_[speciei].dpdv(rho, e, T)
       /sqr(max(rho, this->residualRho_.value()));
}



template<class Thermo>
Foam::scalar
Foam::multicomponentFluidBlastThermo<Thermo>::dpdTi
(
    const label speciei,
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->speciesData_[speciei].dpdT(rho, e, T);
}


template<class Thermo>
Foam::scalar
Foam::multicomponentFluidBlastThermo<Thermo>::mui
(
    const label speciei,
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
    return p;
}


template<class Thermo>
Foam::scalar
Foam::multicomponentFluidBlastThermo<Thermo>::mui
(
    const label speciei,
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    return this->speciesData_[speciei].mu(rho, e, T);
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidBlastThermo<Thermo>::mui
(
    const label speciei,
    const volScalarField& p,
    const volScalarField& T
) const
{
    NotImplemented;
    return p;
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentFluidBlastThermo<Thermo>::mui
(
    const label speciei,
    const volScalarField& rho,
    const volScalarField& e,
    const volScalarField& T
) const
{
    return this->volScalarFieldSpecieProperty
    (
        speciei,
        "mu",
        dimensionSet(1, -1, -1, 0, 0),
        &Thermo::thermoType::mu,
        rho,
        e,
        T
    );
}

// ************************************************************************* //
