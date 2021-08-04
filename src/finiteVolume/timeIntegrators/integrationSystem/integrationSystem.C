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

#include "integrationSystem.H"
#include "timeIntegrator.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::integrationSystem::integrationSystem
(
    const word& name,
    const fvMesh& mesh
)
:
    meshPtr_(&mesh),
    name_(name),
    timeInt_
    (
        meshPtr_->foundObject<timeIntegrator>("globalTimeIntegrator")
      ? &meshPtr_->lookupObjectRef<timeIntegrator>("globalTimeIntegrator")
      : nullptr
    ),
    nSteps_(timeInt_.valid() ? timeInt_->nSteps() : 0),
    oldIs_(timeInt_.valid() ? timeInt_->oldIs() : labelList()),
    nOld_(timeInt_.valid() ? timeInt_->nOld() : 0),
    deltaIs_(timeInt_.valid() ? timeInt_->deltaIs() : labelList()),
    nDelta_(timeInt_.valid() ? timeInt_->nDelta() : 0)
{}

Foam::integrationSystem::integrationSystem()
:
    meshPtr_(nullptr),
    name_(word::null),
    timeInt_(nullptr),
    nSteps_(0),
    oldIs_(0),
    nOld_(0),
    deltaIs_(0),
    nDelta_(0)
{}

void Foam::integrationSystem::set(const word& name, const fvMesh& mesh)
{
    meshPtr_.reset(&mesh);
    name_ = name;
    if (meshPtr_->foundObject<timeIntegrator>("globalTimeIntegrator"))
    {
        timeInt_.reset
        (
            &meshPtr_->lookupObjectRef<timeIntegrator>("globalTimeIntegrator")
        );

        nSteps_ = timeInt_->nSteps();
        oldIs_ = timeInt_->oldIs();
        nOld_ = timeInt_->nOld();
        deltaIs_ = timeInt_->deltaIs();
        nDelta_ = timeInt_->nDelta();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::integrationSystem::~integrationSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::integrationSystem::step() const
{
    return timeInt_.valid() ? timeInt_->step() : 1;
}


Foam::scalarList Foam::integrationSystem::a() const
{
    return timeInt_.valid() ? timeInt_->a() : scalarList(1, 1.0);
}


Foam::scalarList Foam::integrationSystem::b() const
{
    return timeInt_.valid() ? timeInt_->b() : scalarList(1, 1.0);
}


Foam::scalar Foam::integrationSystem::f() const
{
    return timeInt_.valid() ? timeInt_->f() : 1.0;
}


Foam::scalar Foam::integrationSystem::f0() const
{
    return timeInt_.valid() ? timeInt_->f0() : 0.0;
}


bool Foam::integrationSystem::finalStep() const
{
    return timeInt_.valid() ? timeInt_->finalStep() : true;
}


Foam::dimensionedScalar Foam::integrationSystem::time() const
{
    return meshPtr_->time() - meshPtr_->time().deltaT()*(1.0 - f());
}


Foam::dimensionedScalar Foam::integrationSystem::deltaT() const
{
    return meshPtr_->time().deltaT()*f();
}


bool Foam::integrationSystem::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //
