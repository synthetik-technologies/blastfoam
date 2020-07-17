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

#include "RKF45TimeIntegrator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace timeIntegrators
{
    defineTypeNameAndDebug(RKF45, 0);
    addToRunTimeSelectionTable(timeIntegrator, RKF45, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrators::RKF45::RKF45
(
    const fvMesh& mesh
)
:
    timeIntegrator(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrators::RKF45::~RKF45()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeIntegrators::RKF45::setODEFields(integrationSystem& system)
{
    system.setODEFields
    (
        6,
        {true, false, false, false, false, false},
        {true, true, true, true, true, false}
    );
}


void Foam::timeIntegrators::RKF45::integrate()
{
    // Update and store original fields
    Info<< nl << "RKF45: Step 1" << endl;
    this->updateSystems();
    forAll(systems_, i)
    {
        Info<< "Solving " << systems_[i].name() << endl;
        systems_[i].solve(1, {1.0}, {a10});
    }

    // Update and store original fields
    Info<< nl << "RKF45: Step 2" << endl;
    this->updateSystems();
    forAll(systems_, i)
    {
        Info<< "Solving " << systems_[i].name() << endl;
        systems_[i].solve(2, {0.0, 1.0}, {a20, a21});
    }

    // Update and store original fields
    Info<< nl << "RKF45: Step 3" << endl;
    this->updateSystems();
    forAll(systems_, i)
    {
        Info<< "Solving " << systems_[i].name() << endl;
        systems_[i].solve(3, {0.0, 0.0, 1.0}, {a30, a31, a32});
    }

    // Update and store original fields
    Info<< nl << "RKF45: Step 4" << endl;
    this->updateSystems();
    forAll(systems_, i)
    {
        Info<< "Solving " << systems_[i].name() << endl;
        systems_[i].solve
        (
            4,
            {0.0, 0.0, 0.0, 1.0},
            {a40, a41, a42, a43}
        );
    }

    // Update and store original fields
    Info<< nl << "RKF45: Step 5" << endl;
    this->updateSystems();
    forAll(systems_, i)
    {
        Info<< "Solving " << systems_[i].name() << endl;
        systems_[i].solve
        (
            5,
            {0.0, 0.0, 0.0, 0.0, 1.0},
            {a50, a51, a52, a53, a54}
        );
    }

    // Update and store original fields
    Info<< nl << "RKF45: Step 6" << endl;
    this->updateSystems();
    forAll(systems_, i)
    {
        Info<< "Solving " << systems_[i].name() << endl;
        systems_[i].update();
        systems_[i].solve
        (
            6,
            {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {a604, a614, a624, a634, a644, a654}
        );
    }
}
// ************************************************************************* //
