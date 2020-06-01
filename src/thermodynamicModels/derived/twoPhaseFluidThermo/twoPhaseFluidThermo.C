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

#include "twoPhaseFluidThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseFluidThermo::twoPhaseFluidThermo
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
    fluidThermoModel(name, p, rho, e, T, dict, master),
    phases_(dict.lookup("phases")),
    volumeFraction_
    (
        IOobject
        (
            IOobject::groupName("alpha", phases_[0]),
            rho.time().timeName(),
            rho.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        rho.mesh()
    ),
    rho1_
    (
        IOobject
        (
            IOobject::groupName("rho", phases_[0]),
            rho.time().timeName(),
            rho.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        rho.mesh()
    ),
    thermo1_
    (
        fluidThermoModel::New
        (
            phases_[0],
            p_,
            rho1_,
            e_,
            T_,
            dict.subDict(phases_[0]),
            false
        )
    ),
    rho2_
    (
        IOobject
        (
            IOobject::groupName("rho", phases_[1]),
            rho.time().timeName(),
            rho.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        rho.mesh()
    ),
    thermo2_
    (
        fluidThermoModel::New
        (
            phases_[1],
            p_,
            rho2_,
            e_,
            T_,
            dict.subDict(phases_[1]),
            false
        )
    )
{
    //- Force reading of residual values
    thermo1_->read(dict.subDict(phases_[0]));
    thermo2_->read(dict.subDict(phases_[1]));

    // Update total density
    rho_ = volumeFraction_*rho1_ + (1.0 - volumeFraction_)*rho2_;

    this->initialize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseFluidThermo::~twoPhaseFluidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoPhaseFluidThermo::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    thermo1_->solve(stepi, ai, bi);
    thermo2_->solve(stepi, ai, bi);
}


void Foam::twoPhaseFluidThermo::setODEFields
(
    const label nSteps,
    const labelList& oldIs,
    const label& nOld,
    const labelList& deltaIs,
    const label nDelta
)
{

    thermo1_->setODEFields
    (
        nSteps,
        oldIs,
        nOld,
        deltaIs,
        nDelta
    );
    thermo2_->setODEFields
    (
        nSteps,
        oldIs,
        nOld,
        deltaIs,
        nDelta
    );
}


void Foam::twoPhaseFluidThermo::clearODEFields()
{
    thermo1_->clearODEFields();
    thermo2_->clearODEFields();
}


void Foam::twoPhaseFluidThermo::correct()
{
    if (master_)
    {
        T_ = calcT();
        p_ = calcP();
        p_.max(small);
    }

    thermo1_->correct();
    thermo2_->correct();

    if (viscous_)
    {
        // Update transport coefficients
        mu_ =
            volumeFraction_*thermo1_->mu() + (1.0 - volumeFraction_)*thermo2_->mu();
        alpha_ =
            volumeFraction_*thermo1_->alpha()
          + (1.0 - volumeFraction_)*thermo2_->alpha();
    }
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseFluidThermo::speedOfSound() const
{
    volScalarField alphaXi1(volumeFraction_/(thermo1_->Gamma() - 1.0));
    volScalarField alphaXi2((1.0 - volumeFraction_)/(thermo2_->Gamma() - 1.0));

    volScalarField cSqr
    (
        (
            alphaXi1*rho1_*sqr(thermo1_->speedOfSound())
          + alphaXi2*rho2_*sqr(thermo2_->speedOfSound())
        )
        /(alphaXi1 + alphaXi2)/rho_
    );
    cSqr.max(small);
    return sqrt(cSqr);
}


Foam::tmp<Foam::scalarField>
Foam::twoPhaseFluidThermo::speedOfSound(const label patchi) const
{
    scalarField alphaXi1
    (
        volumeFraction_.boundaryField()[patchi]/(thermo1_->Gamma(patchi) - 1.0)
    );
    scalarField alphaXi2
    (
        (1.0 - volumeFraction_.boundaryField()[patchi])
       /(thermo2_->Gamma(patchi) - 1.0)
    );

    scalarField cSqr
    (
        (
            alphaXi1*rho1_.boundaryField()[patchi]
           *sqr(thermo1_->speedOfSound(patchi))
          + alphaXi2*rho2_.boundaryField()[patchi]
           *sqr(thermo2_->speedOfSound(patchi))
        )
        /(rho_.boundaryField()[patchi]*(alphaXi1 + alphaXi2))
    );
    return sqrt(max(cSqr, small));
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidThermo::calcT() const
{
    return
        volumeFraction_*thermo1_->calcT()
      + (1.0 - volumeFraction_)*thermo2_->calcT();
}


Foam::tmp<Foam::scalarField>
Foam::twoPhaseFluidThermo::TRhoE
(
    const scalarField& T,
    const scalarField& e,
    const label patchi
) const
{
    return
        volumeFraction_.boundaryField()[patchi]
       *thermo1_->TRhoE(T, e, patchi)
      + (1.0 - volumeFraction_.boundaryField()[patchi])
       *thermo2_->TRhoE(T, e, patchi);
}


Foam::scalar Foam::twoPhaseFluidThermo::TRhoEi
(
    const scalar& T,
    const scalar& e,
    const label celli
) const
{
    return
        volumeFraction_[celli]*thermo1_->TRhoEi(T, e, celli)
      + (1.0 - volumeFraction_[celli])*thermo2_->TRhoEi(T, e, celli);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidThermo::calcP() const
{
    volScalarField alphaXi1(volumeFraction_/(thermo1_->Gamma() - 1.0));
    volScalarField alphaXi2((1.0 - volumeFraction_)/(thermo2_->Gamma() - 1.0));

    return
        (
            alphaXi1*thermo1_->calcP()
          + alphaXi2*thermo2_->calcP()
        )/(alphaXi1 + alphaXi2);
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseFluidThermo::calce() const
{
    return
        volumeFraction_*thermo1_->calce()
      + (1.0 - volumeFraction_)*thermo2_->calce();
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseFluidThermo::E() const
{
    return
        volumeFraction_*thermo1_->e()
      + (1.0 - volumeFraction_)*thermo2_->e();
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseFluidThermo::e
(
    const volScalarField& rho,
    const volScalarField& e,
    const volScalarField& T
) const
{
    return
        volumeFraction_*thermo1_->e(rho, e, T)
      + (1.0 - volumeFraction_)*thermo2_->e(rho, e, T);
}


Foam::tmp<Foam::scalarField>
Foam::twoPhaseFluidThermo::e
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    return
        volumeFraction_.boundaryField()[patchi]
       *thermo1_->e(rho, e, T, patchi)
      + (1.0 - volumeFraction_.boundaryField()[patchi])
       *thermo2_->e(rho, e, T, patchi);
}



Foam::tmp<Foam::scalarField>
Foam::twoPhaseFluidThermo::e
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const labelList& faceCells
) const
{
    scalarField vf(UIndirectList<scalar>(volumeFraction_(), faceCells)());
    return
        vf*thermo1_->e(rho, e, T, faceCells)
      + (1.0 - vf)*thermo2_->e(rho, e, T, faceCells);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidThermo::Gamma() const
{
    return
        1.0
       /(
            volumeFraction_/thermo1_->Gamma()
          + (1.0 - volumeFraction_)/thermo2_->Gamma()
        );
}


Foam::tmp<Foam::scalarField>
Foam::twoPhaseFluidThermo::Gamma(const label patchi) const
{
    return
        1.0
       /(
            volumeFraction_.boundaryField()[patchi]/thermo1_->Gamma(patchi)
          + (1.0 - volumeFraction_.boundaryField()[patchi])/thermo2_->Gamma(patchi)
        );
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidThermo::ESource() const
{
    return
        volumeFraction_*thermo1_->ESource()
      + (1.0 - volumeFraction_)*thermo2_->ESource();
}



Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidThermo::W() const
{
    return
        volumeFraction_*thermo1_->W()
      + (1.0 - volumeFraction_)*thermo2_->W();
}


Foam::tmp<Foam::scalarField>
Foam::twoPhaseFluidThermo::W(const label patchi) const
{
    return
        1.0
       /(
            volumeFraction_.boundaryField()[patchi]*thermo1_->W(patchi)
          + (1.0 - volumeFraction_.boundaryField()[patchi])*thermo2_->W(patchi)
        );
}


Foam::scalar Foam::twoPhaseFluidThermo::Wi(const label celli) const
{
    return
        volumeFraction_[celli]*thermo1_->Wi(celli)
      + (1.0 - volumeFraction_[celli])*thermo2_->Wi(celli);
}



Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidThermo::Cp() const
{
    return
        volumeFraction_*thermo1_->Cp()
      + (1.0 - volumeFraction_)*thermo2_->Cp();
}



Foam::tmp<Foam::scalarField>
Foam::twoPhaseFluidThermo::Cp(const label patchi) const
{
    return
        1.0
       /(
            volumeFraction_.boundaryField()[patchi]*thermo1_->Cp(patchi)
          + (1.0 - volumeFraction_.boundaryField()[patchi])*thermo2_->Cp(patchi)
        );
}


Foam::tmp<Foam::scalarField>
Foam::twoPhaseFluidThermo::Cp
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    return
        volumeFraction_.boundaryField()[patchi]
       *thermo1_->Cp(rho, e, T, patchi)
      + (1.0 - volumeFraction_.boundaryField()[patchi])
       *thermo2_->Cp(rho, e, T, patchi);
}


Foam::scalar Foam::twoPhaseFluidThermo::Cpi(const label celli) const
{
    return
        volumeFraction_[celli]*thermo1_->Cpi(celli)
      + (1.0 - volumeFraction_[celli])*thermo2_->Cpi(celli);
}



Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidThermo::Cv() const
{
    return
        volumeFraction_*thermo1_->Cv()
      + (1.0 - volumeFraction_)*thermo2_->Cv();
}



Foam::tmp<Foam::scalarField>
Foam::twoPhaseFluidThermo::Cv(const label patchi) const
{
    return
        1.0
       /(
            volumeFraction_.boundaryField()[patchi]*thermo1_->Cv(patchi)
          + (1.0 - volumeFraction_.boundaryField()[patchi])*thermo2_->Cv(patchi)
        );
}


Foam::tmp<Foam::scalarField>
Foam::twoPhaseFluidThermo::Cv
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    return
        volumeFraction_.boundaryField()[patchi]
       *thermo1_->Cv(rho, e, T, patchi)
      + (1.0 - volumeFraction_.boundaryField()[patchi])
       *thermo2_->Cv(rho, e, T, patchi);
}


Foam::scalar Foam::twoPhaseFluidThermo::Cvi(const label celli) const
{
    return
        volumeFraction_[celli]*thermo1_->Cvi(celli)
      + (1.0 - volumeFraction_[celli])*thermo2_->Cvi(celli);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidThermo::CpByCv() const
{
    return
        volumeFraction_*thermo1_->CpByCv()
      + (1.0 - volumeFraction_)*thermo2_->CpByCv();
}


Foam::tmp<Foam::scalarField>
Foam::twoPhaseFluidThermo::CpByCv(const label patchi) const
{
    return
        1.0
       /(
            volumeFraction_.boundaryField()[patchi]*thermo1_->CpByCv(patchi)
          + (1.0 - volumeFraction_.boundaryField()[patchi])*thermo2_->CpByCv(patchi)
        );

}


Foam::tmp<Foam::scalarField>
Foam::twoPhaseFluidThermo::CpByCv
(
    const scalarField& rho,
    const scalarField& e,
    const scalarField& T,
    const label patchi
) const
{
    return
        volumeFraction_.boundaryField()[patchi]
       *thermo1_->CpByCv(rho, e, T, patchi)
      + (1.0 - volumeFraction_.boundaryField()[patchi])
       *thermo2_->CpByCv(rho, e, T, patchi);
}

// ************************************************************************* //
