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
    const scalarField& rhoCells = this->rho_.primitiveField();
    scalarField& heCells = this->heRef();
    scalarField& TCells = this->TRef().primitiveFieldRef();
    scalarField& CpCells = this->CpRef().primitiveFieldRef();
    scalarField& CvCells = this->CvRef().primitiveFieldRef();
    scalarField& kappaCells = this->kappaRef().primitiveFieldRef();
    vectorField* KappaCells =
        !this->isotropic() ? &this->KappaRef() : nullptr;

    const typename Thermo::thermoType1& t1(*this);
    const typename Thermo::thermoType2& t2(*this);
    forAll(rhoCells, celli)
    {
        const scalar x2 = this->cellx(celli);
        const scalar x1 = 1.0 - x2;

        const scalar rhoi = this->rho_[celli];
        scalar& ei = heCells[celli];
        scalar& Ti = TCells[celli];

        if (x2 < this->residualFac_)
        {
            Ti = t1.TRhoE(Ti, rhoi, ei);
            if (Ti < this->TLow_)
            {
                ei = t1.Es(rhoi, ei, this->TLow_);
                Ti = this->TLow_;
            }

            CpCells[celli] = t1.Cp(rhoi, ei, Ti);
            CvCells[celli] = t1.Cv(rhoi, ei, Ti);
            kappaCells[celli] = t1.kappa(rhoi, ei, Ti);

            if (KappaCells)
            {
                (*KappaCells)[celli] = t1.Kappa(rhoi, ei, Ti);
            }
        }
        else if (x1 < this->residualFac_)
        {
            Ti = t2.TRhoE(Ti, rhoi, ei);
            if (Ti < this->TLow_)
            {
                ei = t2.Es(rhoi, ei, this->TLow_);
                Ti = this->TLow_;
            }

            CpCells[celli] = t2.Cp(rhoi, ei, Ti);
            CvCells[celli] = t2.Cv(rhoi, ei, Ti);
            kappaCells[celli] = t2.kappa(rhoi, ei, Ti);

            if (KappaCells)
            {
                (*KappaCells)[celli] = t2.Kappa(rhoi, ei, Ti);
            }
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

            CpCells[celli] =
                t1.Cp(rhoi, ei, Ti)*x1
              + t2.Cp(rhoi, ei, Ti)*x2;
            CvCells[celli] =
                t1.Cv(rhoi, ei, Ti)*x1
              + t2.Cv(rhoi, ei, Ti)*x2;
            kappaCells[celli] =
                t1.kappa(rhoi, ei, Ti)*x1
              + t2.kappa(rhoi, ei, Ti)*x2;

            if (KappaCells)
            {
                (*KappaCells)[celli] =
                    t1.Kappa(rhoi, ei, Ti)*x1
                  + t2.Kappa(rhoi, ei, Ti)*x2;
            }
        }
    }

    const volScalarField::Boundary& rhoBf =
        this->rho_.boundaryFieldRef();

    volScalarField::Boundary& heBf = this->heRef().boundaryFieldRef();
    volScalarField::Boundary& TBf = this->TRef().boundaryFieldRef();

    volScalarField::Boundary& CpBf = this->CpRef().boundaryFieldRef();
    volScalarField::Boundary& CvBf = this->CvRef().boundaryFieldRef();
    volScalarField::Boundary& kappaBf =
        this->kappaRef().boundaryFieldRef();
    volVectorField::Boundary* KappaBf =
        KappaCells
      ? &(this->KappaRef().boundaryFieldRef())
      : nullptr;

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& prho = rhoBf[patchi];

        tmp<scalarField> tpx(this->x(patchi));
        const scalarField px = tpx();

        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& pkappa = kappaBf[patchi];
        fvPatchVectorField* pKappa =
            KappaBf ? &(*KappaBf)[patchi] : nullptr;

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const scalar x2 = px[facei];
                const scalar x1 = 1.0 - x2;

                const scalar rhoi = prho[facei];
                scalar& ei = phe[facei];
                scalar& Ti = pT[facei];

                if (x2 < this->residualFac_)
                {
                    ei = t1.Es(rhoi, ei, Ti);
                    pCp[facei] = t1.Cp(rhoi, ei, Ti);
                    pCv[facei] = t1.Cv(rhoi, ei, Ti);
                    pkappa[facei] = t1.kappa(rhoi, ei, Ti);

                    if (pKappa)
                    {
                        (*pKappa)[facei] = t1.Kappa(rhoi, ei, Ti);
                    }
                }
                else if (x1 < this->residualFac_)
                {
                    ei = t2.Es(rhoi, ei, Ti);
                    pCp[facei] = t2.Cp(rhoi, ei, Ti);
                    pCv[facei] = t2.Cv(rhoi, ei, Ti);
                    pkappa[facei] = t2.kappa(rhoi, ei, Ti);

                    if (pKappa)
                    {
                        (*pKappa)[facei] = t2.Kappa(rhoi, ei, Ti);
                    }
                }
                else
                {
                    ei =
                        t1.Es(rhoi, ei, Ti)*x1
                      + t2.Es(rhoi, ei, Ti)*x2;
                    pCp[facei] =
                        t1.Cp(rhoi, ei, Ti)*x1
                      + t2.Cp(rhoi, ei, Ti)*x2;
                    pCv[facei] =
                        t1.Cv(rhoi, ei, Ti)*x1
                      + t2.Cv(rhoi, ei, Ti)*x2;
                    pkappa[facei] =
                        t1.kappa(rhoi, ei, Ti)*x1
                      + t2.kappa(rhoi, ei, Ti)*x2;

                    if (pKappa)
                    {
                        (*pKappa)[facei] =
                            t1.Kappa(rhoi, ei, Ti)*x1
                          + t2.Kappa(rhoi, ei, Ti)*x2;
                    }
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
                scalar& ei = phe[facei];
                scalar& Ti = pT[facei];

                if (x2 < this->residualFac_)
                {
                    Ti = t1.TRhoE(Ti, rhoi, ei);
                    if (Ti < this->TLow_)
                    {
                        Ti = this->TLow_;
                        ei = t1.Es(rhoi, ei, Ti);
                    }
                    pCp[facei] = t1.Cp(rhoi, ei, Ti);
                    pCv[facei] = t1.Cv(rhoi, ei, Ti);
                    pkappa[facei] = t1.kappa(rhoi, ei, Ti);

                    if (pKappa)
                    {
                        (*pKappa)[facei] = t1.Kappa(rhoi, ei, Ti);
                    }
                }
                else if (x1 < this->residualFac_)
                {
                    Ti = t2.TRhoE(Ti, rhoi, ei);
                    if (Ti < this->TLow_)
                    {
                        Ti = this->TLow_;
                        ei = t2.Es(rhoi, ei, Ti);
                    }
                    pCp[facei] = t2.Cp(rhoi, ei, Ti);
                    pCv[facei] = t2.Cv(rhoi, ei, Ti);
                    pkappa[facei] = t2.kappa(rhoi, ei, Ti);

                    if (pKappa)
                    {
                        (*pKappa)[facei] = t2.Kappa(rhoi, ei, Ti);
                    }
                }
                else
                {
                    Ti =
                        t1.TRhoE(Ti, rhoi, ei)*x1
                      + t2.TRhoE(Ti, rhoi, ei)*x2;
                    if (Ti < this->TLow_)
                    {
                        Ti = this->TLow_;
                        ei =
                            t1.Es(rhoi, ei, Ti)*x1
                          + t2.Es(rhoi, ei, Ti)*x2;
                    }
                    pCp[facei] =
                        t1.Cp(rhoi, ei, Ti)*x1
                      + t2.Cp(rhoi, ei, Ti)*x2;
                    pCv[facei] =
                        t1.Cv(rhoi, ei, Ti)*x1
                      + t2.Cv(rhoi, ei, Ti)*x2;
                    pkappa[facei] =
                        t1.kappa(rhoi, ei, Ti)*x1
                      + t2.kappa(rhoi, ei, Ti)*x2;

                    if (pKappa)
                    {
                        (*pKappa)[facei] =
                            t1.Kappa(rhoi, ei, Ti)*x1
                          + t2.Kappa(rhoi, ei, Ti)*x2;
                    }
                }
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

    dict.readIfPresent("residualActivation", this->residualFac_);
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
void Foam::detonatingSolidBlastThermo<Thermo>::postExplicit()
{
    activation_->postExplicit();
    afterburn_->postExplicit();
}


template<class Thermo>
void Foam::detonatingSolidBlastThermo<Thermo>::postImplicit()
{
    activation_->postImplicit();
    afterburn_->postImplicit();
}


template<class Thermo>
void Foam::detonatingSolidBlastThermo<Thermo>::storeExplicit()
{
    activation_->storeExplicit();
    afterburn_->storeExplicit();
}


template<class Thermo>
void Foam::detonatingSolidBlastThermo<Thermo>::clear()
{
    activation_->clear();
    afterburn_->clear();
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


// ************************************************************************* //
