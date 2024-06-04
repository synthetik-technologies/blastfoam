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
    const scalarField& rhoCells = this->rho_.primitiveFieldRef();
    scalarField& heCells = this->heRef();
    scalarField& TCells = this->TRef().primitiveFieldRef();
    scalarField& CpCells = this->CpRef().primitiveFieldRef();
    scalarField& CvCells = this->CvRef().primitiveFieldRef();
    scalarField& kappaCells = this->kappaRef().primitiveFieldRef();

    forAll(this->rho_, celli)
    {
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

    if (!isotropic())
    {
        vectorField& KappaCells = this->KappaPtr_();

        forAll(KappaCells, celli)
        {
            const scalar rhoi = rhoCells[celli];
            const scalar ei = heCells[celli];
            const scalar Ti = TCells[celli];

            KappaCells[celli] = t.Kappa(rhoi, ei, Ti);
        }

        volVectorField::Boundary& KappaBf =
            this->KappaPtr_().boundaryFieldRef();

        forAll(this->rho_.boundaryField(), patchi)
        {
            const fvPatchScalarField& prho = rhoBf[patchi];
            const fvPatchScalarField& pT = TBf[patchi];
            const fvPatchScalarField& phe = heBf[patchi];

            fvPatchVectorField& pKappa = KappaBf[patchi];
            forAll(prho, facei)
            {
                const scalar rhoi = prho[facei];
                const scalar  ei = phe[facei];
                const scalar Ti = pT[facei];
                pKappa[facei] = t.Kappa(rhoi, ei, Ti);
            }
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
