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

#include "fluidThermoModel.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(fluidThermoModel, 0);
    defineRunTimeSelectionTable(fluidThermoModel, basic);
    defineRunTimeSelectionTable(fluidThermoModel, detonating);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidThermoModel::fluidThermoModel
(
    const word& phaseName,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master
)
:
    basicThermoModel
    (
        phaseName,
        p,
        rho,
        e,
        T,
        dict,
        master
    ),
    regIOobject
    (
        IOobject
        (
            IOobject::groupName("fluidThermo", phaseName),
            e.mesh().time().timeName(),
            e.mesh()
        )
    ),
    mu_
    (
        IOobject
        (
            IOobject::groupName("thermo:mu", phaseName),
            p_.mesh().time().timeName(),
            p_.mesh()
        ),
        p_.mesh(),
        dimensionedScalar(dimDynamicViscosity, 0.0),
        wordList(p_.boundaryField().types().size(), "zeroGradient")
    ),
    viscous_(true)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidThermoModel::~fluidThermoModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluidThermoModel::initialize()
{
    if (!master_)
    {
        return;
    }
    if (max(e_).value() <= 0.0)
    {
        volScalarField e(calce());
        e_ = e;

        //- Force fixed boundaries to be updates
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
    this->correct();

    if (max(mu_).value() < small)
    {
        viscous_ = false;
    }
}


Foam::tmp<Foam::volScalarField> Foam::fluidThermoModel::mu() const
{
    return mu_;
}


Foam::tmp<Foam::scalarField>
Foam::fluidThermoModel::mu(const label patchi) const
{
    return mu_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField> Foam::fluidThermoModel::nu() const
{
    return mu_/rho_;
}


Foam::tmp<Foam::scalarField>
Foam::fluidThermoModel::nu(const label patchi) const
{
    return mu(patchi)/rho_.boundaryField()[patchi];
}


Foam::scalar Foam::fluidThermoModel::nui(const label celli) const
{
    return mu_[celli]/rho_[celli];
}

Foam::tmp<Foam::volScalarField> Foam::fluidThermoModel::alpha() const
{
    return alpha_;
}


Foam::tmp<Foam::scalarField>
Foam::fluidThermoModel::alpha(const label patchi) const
{
    return alpha_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField> Foam::fluidThermoModel::alphaEff
(
    const volScalarField& alphat
) const
{
    return volScalarField::New
    (
        "alphaEff",
        CpByCv()*(alpha_ + alphat)
    );
}


Foam::tmp<Foam::scalarField> Foam::fluidThermoModel::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return CpByCv(patchi)*(alpha(patchi) + alphat);
}


Foam::tmp<Foam::volScalarField> Foam::fluidThermoModel::alphahe() const
{
    return volScalarField::New
    (
        "alphahe",
        CpByCv()*alpha_
    );
}

Foam::tmp<Foam::scalarField>
Foam::fluidThermoModel::alphahe(const label patchi) const
{
    return CpByCv(patchi)*alpha(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::fluidThermoModel::kappa() const
{
    return volScalarField::New
    (
        "kappa",
        Cp()*alpha_
    );
}


Foam::tmp<Foam::scalarField>
Foam::fluidThermoModel::kappa(const label patchi) const
{
    return Cp(patchi)*alpha_.boundaryField()[patchi];
}


Foam::scalar Foam::fluidThermoModel::kappai(const label celli) const
{
    return this->Cpi(celli)*alpha_[celli];
}


Foam::tmp<Foam::volScalarField> Foam::fluidThermoModel::kappaEff
(
    const volScalarField& alphat
) const
{
    return volScalarField::New
    (
        "kappaEff",
        Cp()*(alpha_ + alphat)
    );
}


Foam::tmp<Foam::scalarField>
Foam::fluidThermoModel::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return Cp(patchi)*(alpha_.boundaryField()[patchi] + alphat);
}

// ************************************************************************* //
