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
    forAll(this->rho_, celli)
    {
        const typename Thermo::thermoType& t(this->mixture_[celli]);
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

        scalar Cpi = t.Cp(rhoi, ei, Ti);
        this->CpRef()[celli] = Cpi;
        this->CvRef()[celli] = t.Cv(rhoi, ei, Ti);
        this->alpha()[celli] = t.kappa(rhoi, ei, Ti)/Cpi;
    }

    this->Thermo::baseThermo::T().correctBoundaryConditions();
    this->he().correctBoundaryConditions();

    forAll(this->rho_.boundaryField(), patchi)
    {
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pT =
            this->Thermo::baseThermo::T().boundaryField()[patchi];
        const fvPatchScalarField& phe = this->he().boundaryField()[patchi];

        fvPatchScalarField& pCp = this->CpRef().boundaryFieldRef()[patchi];
        fvPatchScalarField& pCv = this->CvRef().boundaryFieldRef()[patchi];
        fvPatchScalarField& palpha =
            this->alpha().boundaryFieldRef()[patchi];

        forAll(prho, facei)
        {
            const typename Thermo::thermoType& t
            (
                this->mixture_.boundary(patchi, facei)
            );
            const scalar rhoi(prho[facei]);
            const scalar ei(phe[facei]);
            const scalar Ti(pT[facei]);

            const scalar Cpi = t.Cp(rhoi, ei, Ti);
            pCp[facei] = Cpi;
            pCv[facei] = t.Cv(rhoi, ei, Ti);
            palpha[facei] = t.kappa(rhoi, ei, Ti)/Cpi;
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
    volScalarField& rhoRef(this->rho_);
    volScalarField rhoNew
    (
        this->volScalarFieldProperty
        (
            "rho",
            dimDensity,
            &Thermo::thermoType::rho0
        )
    );
    forAll(rhoRef, celli)
    {
        rhoRef[celli] = rhoNew[celli];
    }
    forAll(rhoRef.boundaryField(), patchi)
    {
        forAll(rhoRef.boundaryField()[patchi], facei)
        {
            rhoRef.boundaryFieldRef()[patchi][facei] =
                rhoNew.boundaryField()[patchi][facei];
        }
    }
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::multicomponentSolidBlastThermo<Thermo>::ESource() const
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
Foam::tmp<Foam::volScalarField>
Foam::multicomponentSolidBlastThermo<Thermo>::kappa() const
{
    return this->volScalarFieldProperty
    (
        "kappa",
        dimEnergy/dimTime/dimLength/dimTemperature,
        &Thermo::thermoType::kappa,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class Thermo>
Foam::tmp<Foam::volVectorField>
Foam::multicomponentSolidBlastThermo<Thermo>::Kappa() const
{
    const fvMesh& mesh = this->T_.mesh();

    tmp<volVectorField> tKappa
    (
        volVectorField::New
        (
            "Kappa",
            mesh,
            dimEnergy/dimTime/dimLength/dimTemperature
        )
    );

    volVectorField& Kappa = tKappa.ref();
    vectorField& KappaCells = Kappa.primitiveFieldRef();
    const scalarField& rhoCells = this->rho_;
    const scalarField& eCells = this->e_;
    const scalarField& TCells = this->T_;

    forAll(KappaCells, celli)
    {
        Kappa[celli] =
            this->mixture_[celli].Kappa
            (
                rhoCells[celli],
                eCells[celli],
                TCells[celli]
            );
    }

    volVectorField::Boundary& KappaBf = Kappa.boundaryFieldRef();

    forAll(KappaBf, patchi)
    {
        vectorField& Kappap = KappaBf[patchi];
        const scalarField& pRho = this->rho_.boundaryField()[patchi];
        const scalarField& pe = this->e_.boundaryField()[patchi];
        const scalarField& pT = this->T_.boundaryField()[patchi];

        forAll(Kappap, facei)
        {
            Kappap[facei] =
                this->mixture_.boundary(patchi, facei).Kappa
                (
                    pRho[facei],
                    pe[facei],
                    pT[facei]
                );
        }
    }

    return tKappa;
}


template<class Thermo>
Foam::tmp<Foam::vectorField>
Foam::multicomponentSolidBlastThermo<Thermo>::Kappa(const label patchi) const
{
    const scalarField& pRho = this->rho_.boundaryField()[patchi];
    const scalarField& pe = this->e_.boundaryField()[patchi];
    const scalarField& pT = this->T_.boundaryField()[patchi];
    tmp<vectorField> tKappa(new vectorField(pe.size()));

    vectorField& Kappap = tKappa.ref();

    forAll(Kappap, facei)
    {
        Kappap[facei] =
            this->mixture_.boundary(patchi, facei).Kappa
            (
                pRho[facei],
                pe[facei],
                pT[facei]
            );
    }

    return tKappa;
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
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "mu",
                this->rho_.mesh().time().timeName(),
                this->rho_.mesh()
            ),
            this->rho_.mesh(),
            dimensionedScalar("0", dimensionSet(1, -1, -1, 0, 0), 0.0)
        )
    );
}


// ************************************************************************* //
