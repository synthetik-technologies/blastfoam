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

#include "multicomponentSolidBlastThermo.H"
#include "fvc.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
void Foam::multicomponentSolidBlastThermo<Thermo>::calculate()
{
    this->updateMixture();

    const scalarField& rhoCells = this->rho_.primitiveFieldRef();
    scalarField& heCells = this->heRef();
    scalarField& TCells = this->TRef().primitiveFieldRef();
    scalarField& CpCells = this->CpRef().primitiveFieldRef();
    scalarField& CvCells = this->CvRef().primitiveFieldRef();
    scalarField& kappaCells = this->kappaRef().primitiveFieldRef();

    forAll(this->rho_, celli)
    {
        const typename Thermo::thermoType& t(this->mixture_[celli]);
        const scalar& rhoi(rhoCells[celli]);
        scalar& ei(heCells[celli]);
        scalar& Ti = TCells[celli];

        // Update temperature
        Ti = t.TRhoE(Ti, rhoi, ei);
        if (Ti < this->TLow_)
        {
            ei = t.Es(rhoi, ei, this->TLow_);
            Ti = this->TLow_;
        }

        CpCells[celli] = t.Cp(rhoi, ei, Ti);
        CvCells[celli] = t.Cv(rhoi, ei, Ti);
        kappaCells[celli] = t.kappa(rhoi, ei, Ti);
    }

    const volScalarField::Boundary& rhoBf =
        this->rho_.boundaryFieldRef();

    volScalarField::Boundary& heBf = this->heRef().boundaryFieldRef();
    volScalarField::Boundary& TBf = this->TRef().boundaryFieldRef();

    volScalarField::Boundary& CpBf = this->CpRef().boundaryFieldRef();
    volScalarField::Boundary& CvBf = this->CvRef().boundaryFieldRef();
    volScalarField::Boundary& kappaBf =
        this->kappaRef().boundaryFieldRef();

    forAll(this->rho_.boundaryField(), patchi)
    {
        const fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];

        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(prho, facei)
            {
                const typename Thermo::thermoType& t
                (
                    this->mixture_.boundary(patchi, facei)
                );
                const scalar rhoi(prho[facei]);
                scalar& ei(phe[facei]);
                const scalar Ti(pT[facei]);

                ei = t.Es(rhoi, ei, Ti);

                pCp[facei] = t.Cp(rhoi, ei, Ti);
                pCv[facei] = t.Cv(rhoi, ei, Ti);
                pkappa[facei] = t.kappa(rhoi, ei, Ti);
            }
        }
        else
        {
            forAll(prho, facei)
            {
                const typename Thermo::thermoType& t
                (
                    this->mixture_.boundary(patchi, facei)
                );
                const scalar rhoi(prho[facei]);
                scalar& ei(phe[facei]);
                scalar& Ti(pT[facei]);

                Ti = t.TRhoE(Ti, rhoi, ei);
                if (Ti < this->TLow_)
                {
                    ei = t.Es(rhoi, ei, this->TLow_);
                    Ti = this->TLow_;
                }

                pCp[facei] = t.Cp(rhoi, ei, Ti);
                pCv[facei] = t.Cv(rhoi, ei, Ti);
                pkappa[facei] = t.kappa(rhoi, ei, Ti);
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::multicomponentSolidBlastThermo<Thermo>::multicomponentSolidBlastThermo
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
    this->initializeFields();
}


template<class Thermo>
Foam::multicomponentSolidBlastThermo<Thermo>::multicomponentSolidBlastThermo
(
    const HashPtrTable<Thermo, word, string::hash>& thermoData,
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
    this->initializeFields();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::multicomponentSolidBlastThermo<Thermo>::~multicomponentSolidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::multicomponentSolidBlastThermo<Thermo>::correct()
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
void Foam::multicomponentSolidBlastThermo<Thermo>::updateRho()
{
    this->rho_ == this->volScalarFieldProperty
    (
        "rho",
        dimDensity,
        &Thermo::thermoType::rho0
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentSolidBlastThermo<Thermo>::ESource() const
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
Foam::multicomponentSolidBlastThermo<Thermo>::calce() const
{
    return this->volScalarFieldProperty
    (
        "e",
        dimEnergy/dimMass,
        &Thermo::thermoType::Es,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class Thermo>
Foam::scalar
Foam::multicomponentSolidBlastThermo<Thermo>::cellp(const label celli) const
{
    return this->mixture_[celli].pRhoT
    (
        this->rho_[celli],
        this->e_[celli],
        this->T_[celli]
    );
}


template<class Thermo>
Foam::scalar
Foam::multicomponentSolidBlastThermo<Thermo>::p
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
Foam::multicomponentSolidBlastThermo<Thermo>::p
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
Foam::multicomponentSolidBlastThermo<Thermo>::dpdRho
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
Foam::multicomponentSolidBlastThermo<Thermo>::dpdT
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
Foam::multicomponentSolidBlastThermo<Thermo>::mu
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
Foam::multicomponentSolidBlastThermo<Thermo>::mu
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
Foam::multicomponentSolidBlastThermo<Thermo>::mu
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
Foam::multicomponentSolidBlastThermo<Thermo>::mu
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
