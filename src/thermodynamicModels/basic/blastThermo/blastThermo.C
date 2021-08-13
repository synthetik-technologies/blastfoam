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

#include "blastThermo.H"
#include "basicBlastThermo.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(blastThermo, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blastThermo::blastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    integrationSystem
    (
        IOobject::groupName("blastThermo", phaseName),
        mesh
    ),
    basicThermo::implementation(mesh, dict, phaseName),
    e_
    (
        IOobject
        (
            basicThermo::phasePropertyName
            (
                "e", phaseName
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimMass, 0.0),
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    ),
    rho_
    (
        IOobject
        (
            IOobject::groupName("rho", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimDensity
    ),
    limit_(dict.lookupOrDefault("limit", true))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blastThermo::~blastThermo()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::blastThermo::thermoName() const
{
    return thermo(0).thermoName();
}


const Foam::volScalarField& Foam::blastThermo::he() const
{
    return e_;
}


Foam::volScalarField& Foam::blastThermo::he()
{
    return e_;
}


Foam::tmp<Foam::volScalarField> Foam::blastThermo::rho() const
{
    return rho_;
}


Foam::tmp<Foam::scalarField> Foam::blastThermo::rho(const label patchi) const
{
    return rho_.boundaryField()[patchi];
}


Foam::scalar Foam::blastThermo::rhoi(const label celli) const
{
    return rho_[celli];
}


Foam::volScalarField& Foam::blastThermo::rho()
{
    return rho_;
}


Foam::tmp<Foam::volScalarField> Foam::blastThermo::rho0() const
{
    return rho_.oldTime();
}


bool Foam::blastThermo::limit() const
{
    return limit_;
}


Foam::speciesTable Foam::blastThermo::species() const
{
    return speciesTable();
}


bool Foam::blastThermo::contains(const word&) const
{
    return false;
}


Foam::tmp<Foam::volScalarField> Foam::blastThermo::gamma() const
{
    return this->Cp()/this->Cv();
}


Foam::tmp<Foam::scalarField> Foam::blastThermo::gamma
(
    const scalarField& T,
    const label patchi
) const
{
    return this->Cp(T, patchi)/this->Cv(T, patchi);
}


const Foam::volScalarField& Foam::blastThermo::Y(const word& name) const
{
    NotImplemented;
    return volScalarField::New
    (
        phasePropertyName(name),
        rho_.mesh(),
        dimless
    );
}


Foam::volScalarField& Foam::blastThermo::Y(const word& name)
{
    NotImplemented;
    return rho_;
}


const Foam::volScalarField& Foam::blastThermo::Y(const label i) const
{
    NotImplemented;
    return rho_;
}


Foam::volScalarField& Foam::blastThermo::Y(const label i)
{
    NotImplemented;
    return rho_;
}


void Foam::blastThermo::addDelta
(
    const word&,
    const volScalarField&
)
{}



Foam::tmp<Foam::volScalarField> Foam::blastThermo::kappa() const
{
    return this->Cp()*this->alpha_;
}


Foam::tmp<Foam::scalarField> Foam::blastThermo::kappa
(
    const label patchi
) const
{
    return
        this->Cp(this->T_.boundaryField()[patchi], patchi)
       *this->alpha_.boundaryField()[patchi];
}


Foam::scalar Foam::blastThermo::kappai(const label celli) const
{
    return this->Cpi(celli)*this->alpha_[celli];
}


Foam::tmp<Foam::volScalarField> Foam::blastThermo::alphahe() const
{
    return this->gamma()*this->alpha_;
}


Foam::tmp<Foam::scalarField> Foam::blastThermo::alphahe
(
    const label patchi
) const
{
    return
        this->gamma(this->T_.boundaryField()[patchi], patchi)
       *this->alpha_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField> Foam::blastThermo::kappaEff
(
    const volScalarField& alphat
) const
{
    return this->Cp()*(this->alpha_ + alphat);
}


Foam::tmp<Foam::scalarField> Foam::blastThermo::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        this->Cp(this->T_.boundaryField()[patchi], patchi)
       *(this->alpha_.boundaryField()[patchi] + alphat);
}


Foam::tmp<Foam::volScalarField> Foam::blastThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    return this->gamma()*(this->alpha_ + alphat);
}


Foam::tmp<Foam::scalarField> Foam::blastThermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        this->gamma(this->T_.boundaryField()[patchi], patchi)
       *(this->alpha_.boundaryField()[patchi] + alphat);
}

// ************************************************************************* //
