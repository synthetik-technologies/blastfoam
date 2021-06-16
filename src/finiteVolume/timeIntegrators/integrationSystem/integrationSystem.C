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
    mesh_(mesh),
    name_(name),
    timeInt_
    (
        mesh_.foundObject<timeIntegrator>("globalTimeIntegrator")
      ? &mesh_.lookupObjectRef<timeIntegrator>("globalTimeIntegrator")
      : nullptr
    ),
    nSteps_(timeInt_ ? timeInt_->nSteps() : 0),
    oldIs_(timeInt_ ? timeInt_->oldIs() : labelList()),
    nOld_(timeInt_ ? timeInt_->nOld() : 0),
    deltaIs_(timeInt_ ? timeInt_->deltaIs() : labelList()),
    nDelta_(timeInt_ ? timeInt_->nDelta() : 0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::integrationSystem::~integrationSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::integrationSystem::step() const
{
    return timeInt_ ? timeInt_->step() : 1;
}


Foam::scalarList Foam::integrationSystem::a() const
{
    return timeInt_ ? timeInt_->a() : scalarList(1, 1.0);
}


Foam::scalarList Foam::integrationSystem::b() const
{
    return timeInt_ ? timeInt_->b() : scalarList(1, 1.0);
}


Foam::scalar Foam::integrationSystem::f() const
{
    return timeInt_ ? timeInt_->f() : 1.0;
}


Foam::scalar Foam::integrationSystem::f0() const
{
    return timeInt_ ? timeInt_->f0() : 0.0;
}


bool Foam::integrationSystem::finalStep() const
{
    return timeInt_ ? timeInt_->finalStep() : true;
}


Foam::dimensionedScalar Foam::integrationSystem::time() const
{
    return mesh_.time() - mesh_.time().deltaT()*(1.0 - f());
}


Foam::dimensionedScalar Foam::integrationSystem::deltaT() const
{
    return mesh_.time().deltaT()*f();
}


bool Foam::integrationSystem::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //
