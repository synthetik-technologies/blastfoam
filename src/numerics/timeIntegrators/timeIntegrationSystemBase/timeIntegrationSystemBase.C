/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "timeIntegrationSystemBase.H"
#include "timeIntegrator.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrationSystemBase::timeIntegrationSystemBase
(
    const word& name,
    const objectRegistry& obr
)
:
    name_(name),
    timeInt_
    (
        obr.foundObject<timeIntegrator>(timeIntegrator::typeName)
      ? &obr.lookupObjectRef<timeIntegrator>(timeIntegrator::typeName)
      : nullptr
    )
{}


Foam::timeIntegrationSystemBase::timeIntegrationSystemBase
(
    const word& name,
    const timeIntegrator& integrator
)
:
    name_(name),
    timeInt_(&integrator)
{}



void Foam::timeIntegrationSystemBase::set(const timeIntegrator& integrator)
{
    if (!timeInt_.valid())
    {
        timeInt_.reset(&integrator);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrationSystemBase::~timeIntegrationSystemBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::timeIntegrationSystemBase::step() const
{
    return timeInt_->step();
}


const Foam::scalarList& Foam::timeIntegrationSystemBase::a() const
{
    return timeInt_->a();
}


const Foam::scalarList& Foam::timeIntegrationSystemBase::b() const
{
    return timeInt_->b();
}


Foam::scalar Foam::timeIntegrationSystemBase::f() const
{
    return timeInt_->f();
}


Foam::scalar Foam::timeIntegrationSystemBase::f0() const
{
    return timeInt_->f0();
}


bool Foam::timeIntegrationSystemBase::finalStep() const
{
    return timeInt_->finalStep();
}


Foam::dimensionedScalar Foam::timeIntegrationSystemBase::t() const
{
    return
        timeInt_->time()
      - timeInt_->time().deltaT()*(1.0 - timeInt_->f());
}


Foam::dimensionedScalar Foam::timeIntegrationSystemBase::dt() const
{
    return timeInt_->time().deltaT()*timeInt_->f();
}


bool Foam::timeIntegrationSystemBase::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //
