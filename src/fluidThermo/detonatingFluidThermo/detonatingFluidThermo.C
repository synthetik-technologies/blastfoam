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

#include "detonatingFluidThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class uThermo, class rThermo>
template<class uMethod, class rMethod, class ... Args>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::volScalarFieldProperty
(
    const word& psiName,
    const dimensionSet& psiDim,
    uMethod upsiMethod,
    rMethod rpsiMethod,
    const Args& ... args
) const
{
    tmp<volScalarField> tPsi
    (
        volScalarField::New
        (
            IOobject::groupName(psiName, this->group()),
            this->p_.mesh(),
            psiDim
        )
    );

    volScalarField& psi = tPsi.ref();
    volScalarField x(activation_->lambdaPow());

    forAll(this->p_, celli)
    {
        psi[celli] =
            (this->*rpsiMethod)(args[celli] ...)*x[celli]
          + (this->*upsiMethod)(args[celli] ...)*(1.0 - x[celli]);
    }

    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();

    forAll(psiBf, patchi)
    {
        fvPatchScalarField& pPsi = psiBf[patchi];

        forAll(this->p_.boundaryField()[patchi], facei)
        {
            pPsi[facei] =
                (this->*rpsiMethod)
                (
                    args.boundaryField()[patchi][facei] ...
                )*x.boundaryField()[patchi][facei]
              + (this->*upsiMethod)
                (
                    args.boundaryField()[patchi][facei] ...
                )*(1.0 - x.boundaryField()[patchi][facei]);
        }
    }

    return tPsi;
}


template<class uThermo, class rThermo>
template<class Method, class ... Args>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::singleVolScalarFieldProperty
(
    const word& psiName,
    const dimensionSet& psiDim,
    Method psiMethod,
    const Args& ... args
) const
{
    tmp<volScalarField> tPsi
    (
        volScalarField::New
        (
            IOobject::groupName(psiName, this->group()),
            this->p_.mesh(),
            psiDim
        )
    );

    volScalarField& psi = tPsi.ref();

    forAll(this->p_, celli)
    {
        psi[celli] =
            (this->*psiMethod)(args[celli] ...);
    }

    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();

    forAll(psiBf, patchi)
    {
        fvPatchScalarField& pPsi = psiBf[patchi];

        forAll(this->p_.boundaryField()[patchi], facei)
        {
            pPsi[facei] =
                (this->*psiMethod)
                (
                    args.boundaryField()[patchi][facei] ...
                );
        }
    }

    return tPsi;
}



template<class uThermo, class rThermo>
template<class uMethod, class rMethod, class ... Args>
Foam::tmp<Foam::scalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::cellSetProperty
(
    uMethod upsiMethod,
    rMethod rpsiMethod,
    const labelList& cells,
    const Args& ... args
) const
{
    // Note: Args are fields for the set, not for the mesh as a whole. The
    // cells list is only used to get the mixture.

    tmp<scalarField> tPsi(new scalarField(cells.size()));
    scalarField& psi = tPsi.ref();

    forAll(cells, celli)
    {
        scalar x = activation_->lambdaPowi(celli);
        psi[celli] =
           (this->*rpsiMethod)(args[celli] ...)*x
         + (this->*upsiMethod)(args[celli] ...)*(1.0 - x);
    }

    return tPsi;
}


template<class uThermo, class rThermo>
template<class uMethod, class rMethod, class ... Args>
Foam::tmp<Foam::scalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::patchFieldProperty
(
    uMethod upsiMethod,
    rMethod rpsiMethod,
    const label patchi,
    const Args& ... args
) const
{
    tmp<scalarField> tPsi
    (
        new scalarField(this->p_.boundaryField()[patchi].size())
    );
    scalarField& psi = tPsi.ref();
    scalarField x(activation_->lambdaPow(patchi));

    forAll(this->p_.boundaryField()[patchi], facei)
    {
        psi[facei] =
            (this->*rpsiMethod)(args[facei] ...)*x[facei]
          + (this->*upsiMethod)(args[facei] ...)*(1.0 - x[facei]);
    }

    return tPsi;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class uThermo, class rThermo>
Foam::detonatingFluidThermo<uThermo, rThermo>::detonatingFluidThermo
(
    const word& name,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master
)
:
    fluidThermoModel
    (
        name,
        p,
        rho,
        e,
        T,
        dict,
        master
    ),
    uThermo(dict.subDict("reactants")),
    rThermo(dict.subDict("products")),
    activation_(activationModel::New(rho.mesh(), dict, name)),
    afterburn_(afterburnModel::New(rho.mesh(), dict, name))
{
    //- Initialize the density using the pressure and temperature
    //  This is only done at the first time step (Not on restart)
    if
    (
        uThermo::solid()
     && dict.lookupOrDefault<Switch>("calculateDensity", false)
     && rho_.time().value() == rho_.time().startTime().value()
    )
    {
        volScalarField rhoInit
        (
            volScalarFieldProperty
            (
                IOobject::groupName("rhoInit", name_),
                dimDensity,
                &uThermo::initializeRho,
                &rThermo::initializeRho,
                p_,
                rho_,
                e_,
                T_
            )
        );
        rho_ = rhoInit;
        forAll(rho_.boundaryField(), patchi)
        {
            forAll(rho_.boundaryField()[patchi], facei)
            {
                rho_.boundaryFieldRef()[patchi][facei] =
                    rhoInit.boundaryField()[patchi][facei];
            }
        }
    }

    //- If this is the top level model, initialize the internal energy
    //  if it has not been read
    if (master && max(e_).value() < 0.0)
    {
        volScalarField e(calce());
        e_ = e;
        forAll(e_.boundaryField(), patchi)
        {
            forAll(e_.boundaryField()[patchi], facei)
            {
                e_.boundaryFieldRef()[patchi][facei] =
                    e.boundaryField()[patchi][facei];
            }
        }
        eBoundaryCorrection();
    }
    correct();

    if (max(mu_).value() < small && master)
    {
        viscous_ = false;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class uThermo, class rThermo>
Foam::detonatingFluidThermo<uThermo, rThermo>::~detonatingFluidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class uThermo, class rThermo>
void Foam::detonatingFluidThermo<uThermo, rThermo>::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    activation_->solve(stepi, ai, bi);
    afterburn_->solve(stepi, ai, bi);
}


template<class uThermo, class rThermo>
void Foam::detonatingFluidThermo<uThermo, rThermo>::setODEFields
(
    const label nSteps,
    const labelList& oldIs,
    const label& nOld,
    const labelList& deltaIs,
    const label nDelta
)
{
    activation_->setODEFields
    (
        nSteps,
        oldIs,
        nOld,
        deltaIs,
        nDelta
    );
    afterburn_->setODEFields
    (
        nSteps,
        oldIs,
        nOld,
        deltaIs,
        nDelta
    );
}


template<class uThermo, class rThermo>
void Foam::detonatingFluidThermo<uThermo, rThermo>::clearODEFields()
{
    activation_->clearODEFields();
    afterburn_->clearODEFields();
}


template<class uThermo, class rThermo>
void Foam::detonatingFluidThermo<uThermo, rThermo>::correct()
{

    volScalarField x(activation_->lambdaPow());
    if (master_)
    {
        T_ = calcT();
        p_ = calcP();
        p_.max(small);
    }

    if (viscous_)
    {
        mu_ = volScalarFieldProperty
            (
                IOobject::groupName("mu", name_),
                dimDynamicViscosity,
                &uThermo::mu,
                &rThermo::mu,
                rho_,
                e_,
                T_
            );

        alpha_ = volScalarFieldProperty
            (
                IOobject::groupName("alpha", name_),
                dimensionSet(1, -1, -1, 0, 0),
                &uThermo::alphah,
                &rThermo::alphah,
                rho_,
                e_,
                T_
            );
    }
}


template<class uThermo, class rThermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::speedOfSound() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("speedOfSound", name_),
        dimVelocity,
        &uThermo::speedOfSound,
        &rThermo::speedOfSound,
        p_,
        rho_,
        e_,
        T_
    );
}


template<class uThermo, class rThermo>
Foam::tmp<Foam::scalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::speedOfSound(const label patchi) const
{
    return patchFieldProperty
    (
        &uThermo::speedOfSound,
        &rThermo::speedOfSound,
        patchi,
        p_.boundaryField()[patchi],
        rho_.boundaryField()[patchi],
        e_.boundaryField()[patchi],
        T_.boundaryField()[patchi]
    );
}


template<class uThermo, class rThermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::calcT() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("T", name_),
        dimTemperature,
        &uThermo::TRhoE,
        &rThermo::TRhoE,
        T_,
        rho_,
        e_
    );
}


template<class uThermo, class rThermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::calcP() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("p", name_),
        dimPressure,
        &uThermo::p,
        &rThermo::p,
        rho_,
        e_,
        T_
    );
}


template<class uThermo, class rThermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::calce() const
{
    tmp<volScalarField> eInit
    (
        singleVolScalarFieldProperty
        (
            IOobject::groupName("e", name_),
            dimEnergy/dimMass,
            &uThermo::initializeEnergy,
            p_,
            rho_,
            e_,
            T_
        )
    );

    //- Add detonation energy to initially reacted material
    if (rho_.time().value() == rho_.time().startTime().value())
    {
        eInit.ref() += activation_->e0()*activation_->lambda();
    }

    return eInit;
}


template<class uThermo, class rThermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::E() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("e", name_),
        dimEnergy/dimMass,
        &uThermo::Ea,
        &rThermo::Ea,
        rho_,
        e_,
        T_
    );
}


template<class uThermo, class rThermo>
Foam::tmp<Foam::scalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::E
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &uThermo::Ea,
        &rThermo::Ea,
        patchi,
        rho,
        e,
        T
    );
}


template<class uThermo, class rThermo>
Foam::tmp<Foam::scalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::E
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const labelList& faceCells
) const
{
    return cellSetProperty
    (
        &uThermo::Ea,
        &rThermo::Ea,
        faceCells,
        rho,
        e,
        T
    );
}

template<class uThermo, class rThermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::Gamma() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("Gamma", name_),
        dimless,
        &uThermo::Gamma,
        &rThermo::Gamma,
        rho_,
        e_,
        T_
    );
}


template<class uThermo, class rThermo>
Foam::tmp<Foam::scalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::Gamma(const label patchi) const
{
    return patchFieldProperty
    (
        &uThermo::Gamma,
        &rThermo::Gamma,
        patchi,
        rho_.boundaryField()[patchi],
        e_.boundaryField()[patchi],
        T_.boundaryField()[patchi]
    );
}


template<class uThermo, class rThermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::ESource() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("ESource", name_),
                p_.mesh().time().timeName(),
                p_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            (activation_->ESource() + afterburn_->ESource())*rho_
        )
    );
}


template<class uThermo, class rThermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::W() const
{
    tmp<volScalarField> tPsi
    (
        volScalarField::New
        (
            IOobject::groupName("W", this->group()),
            this->p_.mesh(),
            dimMass/dimMoles
        )
    );

    volScalarField& psi = tPsi.ref();
    volScalarField x(activation_->lambdaPow());

    forAll(this->p_, celli)
    {
        psi[celli] = rThermo::W()*x[celli] + uThermo::W()*(1.0 - x[celli]);
    }

    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();

    forAll(psiBf, patchi)
    {
        fvPatchScalarField& pPsi = psiBf[patchi];

        forAll(this->p_.boundaryField()[patchi], facei)
        {
            pPsi[facei] =
                rThermo::W()*x.boundaryField()[patchi][facei]
              + uThermo::W()*(1.0 - x.boundaryField()[patchi][facei]);
        }
    }

    return tPsi;
}


template<class uThermo, class rThermo>
Foam::scalar
Foam::detonatingFluidThermo<uThermo, rThermo>::Wi(const label celli) const
{
    scalar x = activation_->lambdaPowi(celli);
    return x*uThermo::W() + (1.0 - x)*rThermo::W();
}


template<class uThermo, class rThermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::Cp() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("Cp", name_),
        dimEnergy/dimMass/dimTemperature,
        &uThermo::Cp,
        &rThermo::Cp,
        rho_,
        e_,
        T_
    );
}


template<class uThermo, class rThermo>
Foam::tmp<Foam::scalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::Cp(const label patchi) const
{
    return patchFieldProperty
    (
        &uThermo::Cp,
        &rThermo::Cp,
        patchi,
        rho_.boundaryField()[patchi],
        e_.boundaryField()[patchi],
        T_.boundaryField()[patchi]
    );
}


template<class uThermo, class rThermo>
Foam::scalar
Foam::detonatingFluidThermo<uThermo, rThermo>::Cpi(const label celli) const
{
    scalar x = activation_->lambdaPowi(celli);
    return
        x*uThermo::Cp(rho_[celli], e_[celli], T_[celli])
      + (1.0 - x)*rThermo::Cp(rho_[celli], e_[celli], T_[celli]);
}


template<class uThermo, class rThermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::Cv() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("Cv", name_),
        dimEnergy/dimMass/dimTemperature,
        &uThermo::Cv,
        &rThermo::Cv,
        rho_,
        e_,
        T_
    );
}


template<class uThermo, class rThermo>
Foam::tmp<Foam::scalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::Cv(const label patchi) const
{
    return patchFieldProperty
    (
        &uThermo::Cv,
        &rThermo::Cv,
        patchi,
        rho_.boundaryField()[patchi],
        e_.boundaryField()[patchi],
        T_.boundaryField()[patchi]
    );
}


template<class uThermo, class rThermo>
Foam::tmp<Foam::scalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::Cv
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    return patchFieldProperty
    (
        &uThermo::Cv,
        &rThermo::Cv,
        patchi,
        rho,
        e,
        T
    );
}


template<class uThermo, class rThermo>
Foam::scalar Foam::detonatingFluidThermo<uThermo, rThermo>::Cvi(const label celli) const
{
    scalar x = activation_->lambdaPowi(celli);
    return
        x*uThermo::Cv(rho_[celli], e_[celli], T_[celli])
      + (1.0 - x)*rThermo::Cv(rho_[celli], e_[celli], T_[celli]);
}


template<class uThermo, class rThermo>
Foam::tmp<Foam::volScalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::CpByCv() const
{
    return volScalarFieldProperty
    (
        IOobject::groupName("CpByCv", name_),
        dimless,
        &uThermo::CpByCv,
        &rThermo::CpByCv,
        rho_,
        e_,
        T_
    );
}


template<class uThermo, class rThermo>
Foam::tmp<Foam::scalarField>
Foam::detonatingFluidThermo<uThermo, rThermo>::CpByCv(const label patchi) const
{
    return patchFieldProperty
    (
        &uThermo::CpByCv,
        &rThermo::CpByCv,
        patchi,
        rho_.boundaryField()[patchi],
        e_.boundaryField()[patchi],
        T_.boundaryField()[patchi]
    );
}

// ************************************************************************* //
