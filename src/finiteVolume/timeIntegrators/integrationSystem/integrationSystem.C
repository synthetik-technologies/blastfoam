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
    nSteps_(0),
    timeInt_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::integrationSystem::~integrationSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::integrationSystem::lookupAndInitialize(const word& name)
{
    timeInt_ = &mesh_.lookupObject<timeIntegrator>(name);
    nSteps_ = timeInt_->nSteps();
    timeInt_->setODEFields(*this);
}


void Foam::integrationSystem::setODEFields
(
    const boolList& storeFields,
    const boolList& storeDeltas
)
{
    oldIs_.resize(nSteps_);
    deltaIs_.resize(nSteps_);
    label fi = 0;
    for (label i = 0; i < nSteps_; i++)
    {
        if (storeFields[i])
        {
            oldIs_[i] = fi;
            fi++;
        }
        else
        {
            oldIs_[i] = -1;
        }
    }
    nOld_ = fi;

    fi = 0;
    for (label i = 0; i < nSteps_; i++)
    {
        if (storeDeltas[i])
        {
            deltaIs_[i] = fi;
            fi++;
        }
        else
        {
            deltaIs_[i] = -1;
        }
    }
    nDelta_ = fi;

    this->clearODEFields();
}


Foam::label Foam::integrationSystem::step() const
{
    return timeInt_->step();
}


Foam::scalarList Foam::integrationSystem::a() const
{
    return timeInt_->a();
}


Foam::scalarList Foam::integrationSystem::b() const
{
    return timeInt_->b();
}


Foam::scalar Foam::integrationSystem::f() const
{
    return timeInt_->f();
}


Foam::scalar Foam::integrationSystem::f0() const
{
    return timeInt_->f0();
}


bool Foam::integrationSystem::finalStep() const
{
    return timeInt_->finalStep();
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
