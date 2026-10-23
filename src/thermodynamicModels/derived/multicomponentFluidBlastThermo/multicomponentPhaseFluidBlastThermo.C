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

#include "multicomponentPhaseFluidBlastThermo.H"
#include "fvc.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
void Foam::multicomponentPhaseFluidBlastThermo<Thermo>::calculate
(
    const volScalarField& alpha,
    const volScalarField& he,
    const volScalarField& T,
    volScalarField& alphaCp,
    volScalarField& alphaCv,
    volScalarField& alphamu,
    volScalarField& alphakappa,
    volScalarField& pXiSum,
    volScalarField& XiSum
)
{
    const scalarField& rhoCells = this->rho_.primitiveField();
    const scalarField& alphaCells = alpha.primitiveField();
    const scalarField& heCells = he.primitiveField();
    const scalarField& TCells = T.primitiveField();

    scalarField& alphaCpCells = alphaCp.primitiveFieldRef();
    scalarField& alphaCvCells = alphaCv.primitiveFieldRef();
    scalarField& alphamuCells = alphamu.primitiveFieldRef();
    scalarField& alphakappaCells = alphakappa.primitiveFieldRef();
    scalarField& pXiSumCells = pXiSum.primitiveFieldRef();
    scalarField& XiSumCells = XiSum.primitiveFieldRef();

    forAll(this->rho_, celli)
    {
        if (alphaCells[celli] > 1e-15)
        {
            const typename Thermo::thermoType& t =
                this->cellMixture(celli);

            const scalar alphai = alphaCells[celli];
            const scalar rhoi = rhoCells[celli];
            const scalar ei = heCells[celli];
            const scalar Ti = TCells[celli];

            const scalar Xii = alphai/t.Gamma(rhoi, ei, Ti);

            if (std::isnan(Xii) || std::isinf(Xii) || t.Gamma(rhoi, ei, Ti) < small)
            {
                Pout << "Bad Gamma for multicomponent phase at cell " << celli 
                     << " alpha: " << alphai << " rho: " << rhoi 
                     << " Gamma: " << t.Gamma(rhoi, ei, Ti) << " Xi: " << Xii << endl;
            }

            alphaCpCells[celli] += t.Cp(rhoi, ei, Ti)*alphai*rhoi;
            alphaCvCells[celli] += t.Cv(rhoi, ei, Ti)*alphai*rhoi;
            alphamuCells[celli] += t.mu(rhoi, ei, Ti)*alphai*rhoi;
            alphakappaCells[celli] += t.kappa(rhoi, ei, Ti)*alphai*rhoi;
            pXiSumCells[celli] += t.p(rhoi, ei, Ti)*Xii;
            XiSumCells[celli] += Xii;
        }
    }

    forAll(alpha.boundaryField(), patchi)
    {
        const fvPatchScalarField& palpha = alpha.boundaryField()[patchi];
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& phe = he.boundaryField()[patchi];
        const fvPatchScalarField& pT = T.boundaryField()[patchi];

        fvPatchScalarField& palphaCp = alphaCp.boundaryFieldRef()[patchi];
        fvPatchScalarField& palphaCv = alphaCv.boundaryFieldRef()[patchi];
        fvPatchScalarField& palphamu = alphamu.boundaryFieldRef()[patchi];
        fvPatchScalarField& palphakappa = alphakappa.boundaryFieldRef()[patchi];
        fvPatchScalarField& ppXiSum = pXiSum.boundaryFieldRef()[patchi];
        fvPatchScalarField& pXiSum = XiSum.boundaryFieldRef()[patchi];

        forAll(palpha, facei)
        {
            if (palpha[facei] > 1e-15)
            {
                const typename Thermo::thermoType& t =
                    this->patchFaceMixture(patchi, facei);

                const scalar alphai = palpha[facei];
                const scalar rhoi = prho[facei];
                const scalar ei = phe[facei];
                const scalar Ti = pT[facei];

                const scalar Xii = alphai/t.Gamma(rhoi, ei, Ti);

                palphaCp[facei] += t.Cp(rhoi, ei, Ti)*alphai*rhoi;
                palphaCv[facei] += t.Cv(rhoi, ei, Ti)*alphai*rhoi;
                palphamu[facei] += t.mu(rhoi, ei, Ti)*alphai*rhoi;
                palphakappa[facei] += t.kappa(rhoi, ei, Ti)*alphai*rhoi;
                ppXiSum[facei] += t.p(rhoi, ei, Ti)*Xii;
                pXiSum[facei] += Xii;
            }
        }
    }
}


template<class Thermo>
void Foam::multicomponentPhaseFluidBlastThermo<Thermo>::calculateSpeedOfSound
(
    const volScalarField& alpha,
    volScalarField& cSqrRhoXiSum
)
{
    const scalarField& alphaCells = alpha.primitiveField();
    scalarField& cSqrRhoXiSumCells = cSqrRhoXiSum.primitiveFieldRef();

    forAll(this->rho_, celli)
    {
        const scalar alphai = alphaCells[celli];

        if (alphai > 1e-15)
        {
            const typename Thermo::thermoType& t =
                this->cellMixture(celli);

            const scalar rhoi = this->rho_[celli];
            const scalar ei = this->he()[celli];
            const scalar Ti = this->T()[celli];
            const scalar pi = t.p(rhoi, ei, Ti);

            cSqrRhoXiSumCells[celli] +=
                max(t.cSqr(pi, rhoi, ei, Ti), small)
               *rhoi*alphai/t.Gamma(rhoi, ei, Ti);
        }
    }

    forAll(alpha.boundaryField(), patchi)
    {
        const fvPatchScalarField& palpha = alpha.boundaryField()[patchi];
        fvPatchScalarField& pcSqrRhoXiSum =
            cSqrRhoXiSum.boundaryFieldRef()[patchi];

        forAll(palpha, facei)
        {
            const scalar alphai = palpha[facei];

            if (alphai > 1e-15)
            {
                const typename Thermo::thermoType& t =
                    this->patchFaceMixture(patchi, facei);

                const scalar rhoi = this->rho_.boundaryField()[patchi][facei];
                const scalar ei = this->he().boundaryField()[patchi][facei];
                const scalar Ti = this->T().boundaryField()[patchi][facei];
                const scalar pi = t.p(rhoi, ei, Ti);

                pcSqrRhoXiSum[facei] +=
                    max(t.cSqr(pi, rhoi, ei, Ti), small)
               *rhoi*alphai/t.Gamma(rhoi, ei, Ti);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::multicomponentPhaseFluidBlastThermo<Thermo>::multicomponentPhaseFluidBlastThermo
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
{}


template<class Thermo>
Foam::multicomponentPhaseFluidBlastThermo<Thermo>::multicomponentPhaseFluidBlastThermo
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
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::multicomponentPhaseFluidBlastThermo<Thermo>::~multicomponentPhaseFluidBlastThermo()
{}


// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentPhaseFluidBlastThermo<Thermo>::ESource() const
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
Foam::multicomponentPhaseFluidBlastThermo<Thermo>::ESource
(
    const volScalarField& alpha
) const
{
    return ESource();
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentPhaseFluidBlastThermo<Thermo>::initESource() const
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
Foam::multicomponentPhaseFluidBlastThermo<Thermo>::initESource
(
    const volScalarField& alpha
) const
{
    return initESource();
}


template<class Thermo>
bool Foam::multicomponentPhaseFluidBlastThermo<Thermo>::inviscid() const
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
void Foam::multicomponentPhaseFluidBlastThermo<Thermo>::updateRho
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
Foam::multicomponentPhaseFluidBlastThermo<Thermo>::calce
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
Foam::multicomponentPhaseFluidBlastThermo<Thermo>::calcCelle
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
Foam::scalar Foam::multicomponentPhaseFluidBlastThermo<Thermo>::cellpRhoT
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
Foam::scalar Foam::multicomponentPhaseFluidBlastThermo<Thermo>::patchFacepRhoT
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
Foam::multicomponentPhaseFluidBlastThermo<Thermo>::Gamma() const
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
Foam::multicomponentPhaseFluidBlastThermo<Thermo>::cellGamma(const label celli) const
{
    return this->cellMixture(celli).Gamma
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::scalar Foam::multicomponentPhaseFluidBlastThermo<Thermo>::patchFaceGamma
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
Foam::tmp<Foam::volScalarField>
Foam::multicomponentPhaseFluidBlastThermo<Thermo>::calcSpeedOfSound() const
{
    tmp<volScalarField> tcSqr
    (
        Thermo::volScalarFieldProperty
        (
            "cSqr",
            sqr(dimVelocity),
            &Thermo::thermoType::cSqr,
            this->p(),
            this->rho_,
            this->e_,
            this->T_
        )
    );
    tcSqr.ref().max(small);
    return sqrt(tcSqr);
}


template<class Thermo>
Foam::scalar Foam::multicomponentPhaseFluidBlastThermo<Thermo>::cellSpeedOfSound
(
    const scalar p,
    const label celli
) const
{
    return sqrt
    (
        max
        (
            this->cellMixture(celli).cSqr
            (
                p,
                this->rho_[celli],
                this->e_[celli],
                this->T_[celli]
            ),
            small
        )
    );
}


template<class Thermo>
Foam::scalar
Foam::multicomponentPhaseFluidBlastThermo<Thermo>::celldpdRho(const label celli) const
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
Foam::multicomponentPhaseFluidBlastThermo<Thermo>::celldpde(const label celli) const
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
Foam::multicomponentPhaseFluidBlastThermo<Thermo>::celldpdT(const label celli) const
{
    return this->cellMixture(celli).dpdT
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


// ************************************************************************* //
