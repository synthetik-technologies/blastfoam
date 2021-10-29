/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "AnisothermalImmersedBoundaryObject.H"
#include "fvm.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ImmersedType>
Foam::AnisothermalImmersedBoundaryObject<ImmersedType>::
AnisothermalImmersedBoundaryObject
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    ImmersedType(mesh, dict),
    thermoDict_
    (
        IOobject
        (
            "thermophysicalProperties",
            mesh.time().constant(),
            *(this->immersedFvMesh()),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    thermo_(),
    kappaPtr_(),
    CpPtr_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ImmersedType>
Foam::AnisothermalImmersedBoundaryObject<ImmersedType>::
~AnisothermalImmersedBoundaryObject()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ImmersedType>
void Foam::AnisothermalImmersedBoundaryObject<ImmersedType>::initialize()
{
    ImmersedType::initialize();
    thermo_.set
    (
        solidThermoModel::New
        (
            word::null,
            *(this->immersedFvMesh()),
            thermoDict_,
            true
        ).ptr()
    );
    kappaPtr_.set(new volScalarField(thermo_->kappa()));
    CpPtr_.set(new volScalarField(thermo_->Cp()));
}

template<class ImmersedType>
void Foam::AnisothermalImmersedBoundaryObject<ImmersedType>::solve
(
    const bool restart
)
{
    const volScalarField& rho(thermo_->rho());
    kappaPtr_() = thermo_->kappa();
    CpPtr_() = thermo_->Cp();
    Foam::solve
    (
        fvm::ddt(rho, CpPtr_(), thermo_->T())
     ==
        fvm::laplacian(kappaPtr_(), thermo_->T())
    );

    ImmersedType::update();
}


template<class ImmersedType>
void
Foam::AnisothermalImmersedBoundaryObject<ImmersedType>::addThermalForcing
(
    const word& name,
    volScalarField& Q,
    const volScalarField& alphaRho,
    const volScalarField& alphaRhoEOld,
    const volScalarField& RHS,
    const dimensionedScalar& dT
) const
{
    this->scalarBoundaries_[name].addForcing
    (
        Q,
        alphaRho,
        alphaRhoEOld,
        RHS,
        dT.value()
    );
}


template<class ImmersedType>
void Foam::AnisothermalImmersedBoundaryObject<ImmersedType>::setValues()
{
    forAll(energyBoundaries_, i)
    {
        energyBoundaries_[i].setValues();
    }
}


template<class ImmersedType>
void Foam::AnisothermalImmersedBoundaryObject<ImmersedType>::status() const
{
    ImmersedType::status();
    const volScalarField& T(thermo_->T());

    Info<< "    Temperature (min, max, avg): "
        << min(T).value() << ", "
        << max(T).value() << ", "
        << T.weightedAverage(this->immersedFvMesh()->V()).value()
        << endl;
}


template<class ImmersedType>
void Foam::AnisothermalImmersedBoundaryObject<ImmersedType>::movePoints()
{
    ImmersedType::movePoints();
}


template<class ImmersedType>
bool Foam::AnisothermalImmersedBoundaryObject<ImmersedType>::read
(
    const dictionary& dict
)
{
    return ImmersedType::read(dict);
}


template<class ImmersedType>
void Foam::AnisothermalImmersedBoundaryObject<ImmersedType>::write
(
    Ostream& os
) const
{
    ImmersedType::write(os);
}


template<class ImmersedType>
void Foam::AnisothermalImmersedBoundaryObject<ImmersedType>::write
(
    dictionary& dict
) const
{
    ImmersedType::write(dict);
}

// ************************************************************************* //
