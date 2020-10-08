/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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

Foam::scalar Foam::timeIntegrator::f() const
{
    return sum(bs_[stepi_]);
}

Foam::scalar Foam::timeIntegrator::f0() const
{
    if (stepi_ == 0)
    {
        return 0.0;
    }

    scalarList t0s(stepi_, 0.0);
    scalarList dts(stepi_, 0.0);
    dts[0] = bs_[0][0];

    for (label i = 1; i < stepi_; i++)
    {
        forAll(as_[i], j)
        {
            t0s[i] += as_[i][j]*t0s[j];
            dts[i] += bs_[i][j];
        }
        t0s[i] += dts[i-1];
    }
    Info<<t0s<<endl;
    return t0s.last();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrator::timeIntegrator(const fvMesh& mesh)
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
    stepi_(0)
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
        6,
        saveOlds,
        saveDeltas
    );
}


void Foam::timeIntegrator::integrate()
{
    // Update and store original fields
    for (stepi_ = 0; stepi_ < as_.size(); stepi_++)
    {
        Info<< nl << this->type() << ": step " << stepi_ << endl;
        this->updateAll();
        forAll(systems_, i)
        {
            Info<< "Solving " << systems_[i].name() << endl;
            systems_[i].solve(stepi_+1, as_[stepi_], bs_[stepi_]);
        }
    }

    this->postUpdateAll();
}


// ************************************************************************* //
