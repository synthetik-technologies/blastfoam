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

#include "basicSolidBlastThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
void Foam::basicSolidBlastThermo<Thermo>::calculate()
{
    const typename Thermo::thermoType& t(*this);
    forAll(this->rho_, celli)
    {
        const scalar& rhoi(this->rho_[celli]);
        scalar& ei(this->heRef()[celli]);
        scalar& Ti = this->TRef()[celli];

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
        this->alphaRef()[celli] = t.kappa(rhoi, ei, Ti)/Cpi;
    }

    this->TRef().correctBoundaryConditions();
    this->heRef().correctBoundaryConditions();

    forAll(this->rho_.boundaryField(), patchi)
    {
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pT =
            this->TRef().boundaryField()[patchi];
        const fvPatchScalarField& phe =
            this->heRef().boundaryField()[patchi];

        fvPatchScalarField& pCp = this->CpRef().boundaryFieldRef()[patchi];
        fvPatchScalarField& pCv = this->CvRef().boundaryFieldRef()[patchi];
        fvPatchScalarField& palpha =
            this->alphaRef().boundaryFieldRef()[patchi];

        forAll(prho, facei)
        {
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
Foam::basicSolidBlastThermo<Thermo>::basicSolidBlastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
:
    Thermo(mesh, dict, phaseName, masterName)
{
    this->initializeFields();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::basicSolidBlastThermo<Thermo>::~basicSolidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::basicSolidBlastThermo<Thermo>::correct()
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
void Foam::basicSolidBlastThermo<Thermo>::updateRho()
{
    this->rho_ == Thermo::volScalarFieldProperty
    (
        "rho",
        dimDensity,
        &Thermo::thermoType::rho0
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::basicSolidBlastThermo<Thermo>::ESource() const
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
Foam::basicSolidBlastThermo<Thermo>::kappa() const
{
    return Thermo::volScalarFieldProperty
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
Foam::basicSolidBlastThermo<Thermo>::Kappa() const
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

    const scalarField& rhoCells = this->rho_;
    const scalarField& eCells = this->e_;
    const scalarField& TCells = this->T_;

    volVectorField& Kappa = tKappa.ref();
    vectorField& KappaCells = Kappa.primitiveFieldRef();

    forAll(KappaCells, celli)
    {
        Kappa[celli] =
            Thermo::thermoType::Kappa(rhoCells[celli], eCells[celli], TCells[celli]);
    }

    volVectorField::Boundary& KappaBf = Kappa.boundaryFieldRef();

    forAll(KappaBf, patchi)
    {
        const scalarField& pRho = this->rho_.boundaryField()[patchi];
        const scalarField& pe = this->e_.boundaryField()[patchi];
        const scalarField& pT = this->T_.boundaryField()[patchi];

        vectorField& Kappap = KappaBf[patchi];

        forAll(Kappap, facei)
        {
            Kappap[facei] =
                Thermo::thermoType::Kappa(pRho[facei], pe[facei], pT[facei]);
        }
    }

    return tKappa;
}


template<class Thermo>
Foam::tmp<Foam::vectorField>
Foam::basicSolidBlastThermo<Thermo>::Kappa
(
    const label patchi
) const
{
    const scalarField& pRho = this->rho_.boundaryField()[patchi];
    const scalarField& pe = this->e_.boundaryField()[patchi];
    const scalarField& pT = this->T_.boundaryField()[patchi];

    tmp<vectorField> tKappa(new vectorField(pT.size()));
    vectorField& Kappap = tKappa.ref();

    forAll(pT, facei)
    {
        Kappap[facei] =
            Thermo::thermoType::Kappa(pRho[facei], pe[facei], pT[facei]);
    }

    return tKappa;
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::basicSolidBlastThermo<Thermo>::calce() const
{
    return Thermo::volScalarFieldProperty
    (
        "e",
        dimEnergy/dimMass,
        &Thermo::Es,
        this->rho_,
        this->e_,
        this->T_
    );
}

// ************************************************************************* //
