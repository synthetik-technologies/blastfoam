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

#include "detonatingSolidBlastThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo>
void Foam::detonatingSolidBlastThermo<Thermo>::calculate()
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

            const scalar Cpi = t1.Cp(rhoi, ei, Ti);
            this->CpRef()[celli] = Cpi;
            this->CvRef()[celli] = t1.Cv(rhoi, ei, Ti);
            this->alphaRef()[celli] = t1.kappa(rhoi, ei, Ti)/Cpi;
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

            const scalar Cpi = t2.Cp(rhoi, ei, Ti);
            this->CpRef()[celli] = Cpi;
            this->CvRef()[celli] = t2.Cv(rhoi, ei, Ti);
            this->alphaRef()[celli] = t2.kappa(rhoi, ei, Ti)/Cpi;
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

            this->CpRef()[celli] =
                t1.Cp(rhoi, ei, Ti)*x1
              + t2.Cp(rhoi, ei, Ti)*x2;;
            this->CvRef()[celli] =
                t1.Cv(rhoi, ei, Ti)*x1
              + t2.Cv(rhoi, ei, Ti)*x2;
            this->alphaRef()[celli] =
                t1.kappa(rhoi, ei, Ti)/t1.Cp(rhoi, ei, Ti)*x1
              + t2.kappa(rhoi, ei, Ti)/t2.Cp(rhoi, ei, Ti)*x2;
        }
    }

    this->TRef().correctBoundaryConditions();
    this->heRef().correctBoundaryConditions();

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& prho = this->rho_.boundaryField()[patchi];
        const fvPatchScalarField& pT =
            this->TRef().boundaryField()[patchi];
        const fvPatchScalarField& phe =
            this->heRef().boundaryField()[patchi];
        const scalarField px(this->x(patchi));

        fvPatchScalarField& pCp = this->CpRef().boundaryFieldRef()[patchi];
        fvPatchScalarField& pCv = this->CvRef().boundaryFieldRef()[patchi];
        fvPatchScalarField& palpha =
            this->alphaRef().boundaryFieldRef()[patchi];

        forAll(pT, facei)
        {
            const scalar x2 = px[facei];
            const scalar x1 = 1.0 - x2;
            const scalar rhoi(prho[facei]);
            const scalar ei(phe[facei]);
            const scalar Ti(pT[facei]);

            if (x2 < this->residualActivation_)
            {
                pCp[facei] = t1.Cp(rhoi, ei, Ti);
                pCv[facei] = t1.Cv(rhoi, ei, Ti);
                palpha[facei] = t1.kappa(rhoi, ei, Ti)/pCp[facei];
            }
            else if (x1 < this->residualActivation_)
            {
                pCp[facei] = t2.Cp(rhoi, ei, Ti);
                pCv[facei] = t2.Cv(rhoi, ei, Ti);
                palpha[facei] = t2.kappa(rhoi, ei, Ti)/pCp[facei];
            }
            else
            {
                pCp[facei] =
                    t1.Cp(rhoi, ei, Ti)*x1
                  + t2.Cp(rhoi, ei, Ti)*x2;
                pCv[facei] =
                    t1.Cv(rhoi, ei, Ti)*x1
                  + t2.Cv(rhoi, ei, Ti)*x2;
                palpha[facei] =
                    t1.kappa(rhoi, ei, Ti)/pCp[facei]*x1
                  + t2.kappa(rhoi, ei, Ti)/pCp[facei]*x2;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::detonatingSolidBlastThermo<Thermo>::detonatingSolidBlastThermo
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
    this->initializeFields();
}


template<class Thermo>
void Foam::detonatingSolidBlastThermo<Thermo>::initializeModels()
{
    activation_->initializeModels();
    afterburn_->initializeModels();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo>
Foam::detonatingSolidBlastThermo<Thermo>::~detonatingSolidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::detonatingSolidBlastThermo<Thermo>::correct()
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
void Foam::detonatingSolidBlastThermo<Thermo>::update()
{
    activation_->update();
    afterburn_->update();
}


template<class Thermo>
void Foam::detonatingSolidBlastThermo<Thermo>::solve()
{
    activation_->solve();
    afterburn_->solve();
}


template<class Thermo>
void Foam::detonatingSolidBlastThermo<Thermo>::postUpdate()
{
    activation_->postUpdate();
    afterburn_->postUpdate();
}



template<class Thermo>
void Foam::detonatingSolidBlastThermo<Thermo>::updateRho()
{
    this->rho_ == Thermo::blendedVolScalarFieldProperty
    (
        "rho",
        dimDensity,
        &Thermo::thermoType1::rho0,
        &Thermo::thermoType2::rho0
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingSolidBlastThermo<Thermo>::calce() const
{
    return volScalarField::New
    (
        "eInit",
        Thermo::blendedVolScalarFieldProperty
        (
            "e",
            dimEnergy/dimMass,
            &Thermo::thermoType1::Es,
            &Thermo::thermoType2::Es,
            this->rho_,
            this->e_,
            this->T_
        ) + activation_->initESource()
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingSolidBlastThermo<Thermo>::ESource() const
{
    return volScalarField::New
    (
        "ESource",
        (activation_->ESource() + afterburn_->ESource())*this->rho_
    );
}


template<class Thermo>
Foam::tmp<Foam::volScalarField> Foam::detonatingSolidBlastThermo<Thermo>::kappa() const
{
    return Thermo::blendedVolScalarFieldProperty
    (
        "kappa",
        dimEnergy/dimTime/dimLength/dimTemperature,
        &Thermo::thermoType1::kappa,
        &Thermo::thermoType2::kappa,
        this->rho_,
        this->e_,
        this->T_
    );
}


template<class Thermo>
Foam::tmp<Foam::volVectorField>
Foam::detonatingSolidBlastThermo<Thermo>::Kappa() const
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
        scalar x = cellx(celli);
        if (x < this->residualActivation_)
        {
            Kappa[celli] =
                Thermo::thermoType1::Kappa
                (
                    rhoCells[celli],
                    eCells[celli],
                    TCells[celli]
                );
        }
        else if ((1.0 - x) < this->residualActivation_)
        {
            Kappa[celli] =
                Thermo::thermoType2::Kappa
                (
                    rhoCells[celli],
                    eCells[celli],
                    TCells[celli]
                );
        }
        else
        {
            Kappa[celli] =
                Thermo::thermoType2::Kappa
                (
                    rhoCells[celli],
                    eCells[celli],
                    TCells[celli]
                )*x
              + Thermo::thermoType1::Kappa
                (
                    rhoCells[celli],
                    eCells[celli],
                    TCells[celli]
                )*(1.0 - x);
        }
    }

    volVectorField::Boundary& KappaBf = Kappa.boundaryFieldRef();

    forAll(KappaBf, patchi)
    {
        vectorField& Kappap = KappaBf[patchi];
        const scalarField& pRho = this->rho_.boundaryField()[patchi];
        const scalarField& pe = this->e_.boundaryField()[patchi];
        const scalarField& pT = this->T_.boundaryField()[patchi];
        tmp<scalarField> xp(this->x(patchi));

        forAll(Kappap, facei)
        {
            const scalar& x = xp()[facei];
            if (x < this->residualActivation_)
            {
                Kappap[facei] =
                    Thermo::thermoType1::Kappa
                    (
                        pRho[facei],
                        pe[facei],
                        pT[facei]
                    );
            }
            else if ((1.0 - x) < this->residualActivation_)
            {
                Kappap[facei] =
                    Thermo::thermoType2::Kappa
                    (
                        pRho[facei],
                        pe[facei],
                        pT[facei]
                    );
            }
            else
            {
                Kappap[facei] =
                    Thermo::thermoType2::Kappa
                    (
                        pRho[facei],
                        pe[facei],
                        pT[facei]
                    )*x
                  + Thermo::thermoType1::Kappa
                    (
                        pRho[facei],
                        pe[facei],
                        pT[facei]
                    )*(1.0 - x);
            }
        }
    }

    return tKappa;
}


template<class Thermo>
Foam::tmp<Foam::vectorField>
Foam::detonatingSolidBlastThermo<Thermo>::Kappa(const label patchi) const
{
    const scalarField& pRho = this->rho_.boundaryField()[patchi];
    const scalarField& pe = this->e_.boundaryField()[patchi];
    const scalarField& pT = this->T_.boundaryField()[patchi];
    tmp<vectorField> tKappa(new vectorField(pe.size()));

    vectorField& Kappap = tKappa.ref();
    tmp<scalarField> xp(this->x(patchi));

    forAll(pe, facei)
    {
        const scalar& x = xp()[facei];
        if (x < this->residualActivation_)
        {
            Kappap[facei] =
                Thermo::thermoType1::Kappa
                (
                    pRho[facei],
                    pe[facei],
                    pT[facei]
                );
        }
        else if ((1.0 - x) < this->residualActivation_)
        {
            Kappap[facei] =
                Thermo::thermoType2::Kappa
                (
                    pRho[facei],
                    pe[facei],
                    pT[facei]
                );
        }
        else
        {
            Kappap[facei] =
                Thermo::thermoType2::Kappa
                (
                    pRho[facei],
                    pe[facei],
                    pT[facei]
                )*x
                + Thermo::thermoType1::Kappa
                (
                    pRho[facei],
                    pe[facei],
                    pT[facei]
                )*(1.0 - x);
        }
    }

    return tKappa;
}


// ************************************************************************* //
