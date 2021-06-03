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

#include "timeIntegrator.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(timeIntegrator, 0);
    defineRunTimeSelectionTable(timeIntegrator, dictionary);
}


// * * * * * * * * * * * * * * Protected Functions * * * * * * * * * * * * * //

void Foam::timeIntegrator::updateAll()
{
    forAll(systems_, i)
    {
        systems_[i].update();
    }
}


void Foam::timeIntegrator::postUpdateAll()
{
    forAll(systems_, i)
    {
        systems_[i].postUpdate();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrator::timeIntegrator(const fvMesh& mesh, const label)
:
    regIOobject
    (
        IOobject
        (
            "globalTimeIntegrator",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    stepi_(0),
    f_(0),
    f0_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrator::~timeIntegrator()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

void Foam::timeIntegrator::addSystem(integrationSystem& system)
{

    label oldSize = systems_.size();

    setODEFields(system);
    systems_.resize(oldSize + 1);
    systems_.set(oldSize, &system);

    if (f_.size() == 0)
    {
        f_.resize(nSteps());
        f0_.resize(nSteps());
        f0_[0] = 0.0;
        f_[0] = sum(bs_[0]);
        for (label stepi = 2; stepi <= nSteps(); stepi++)
        {
            f_[stepi-1] = sum(bs_[stepi-1]);
            scalarList ts(stepi+1, 0.0);
            scalarList dts(stepi, 0.0);
            forAll(dts, i)
            {
                dts[i] = sum(bs_[i]);
            }
            ts[1] = dts[0];

            for (label i = 1; i < stepi; i++)
            {
                for (label j = 0; j < as_[i].size(); j++)
                {
                    ts[i+1] += as_[i][j]*ts[j];
                }
                ts[i+1] += dts[i];
            }
            f0_[stepi-1] = ts.last() - dts.last();
        }
    }
}


void Foam::timeIntegrator::setODEFields(integrationSystem& system) const
{
    //- Determine if what fields need to be saved
    boolList saveOlds(as_.size(), false);
    boolList saveDeltas(bs_.size(), false);

    forAll(as_, i)
    {
        for (label j = 0; j < as_[i].size() - 1; j++)
        {
            saveOlds[j] = saveOlds[j] || mag(as_[i][j]) > small;
            saveDeltas[j] = saveDeltas[j] || mag(bs_[i][j]) > small;
        }
    }

    system.setODEFields
    (
        saveOlds,
        saveDeltas
    );
}


void Foam::timeIntegrator::integrate()
{
    // Update and store original fields
    for (stepi_ = 1; stepi_ <= as_.size(); stepi_++)
    {
        Info<< nl << this->type() << ": step " << stepi_ << endl;
        this->updateAll();
        forAll(systems_, i)
        {
            Info<< "Solving " << systems_[i].name() << endl;
            systems_[i].solve();
        }
    }

    this->postUpdateAll();
    stepi_ = 0;
}


Foam::scalar Foam::timeIntegrator::f() const
{
    return f_[stepi_-1];
}

Foam::scalar Foam::timeIntegrator::f0() const
{
    return f0_[stepi_ - 1];
}


bool Foam::timeIntegrator::finalStep() const
{
    return stepi_ == as_.size();
}



// ************************************************************************* //
