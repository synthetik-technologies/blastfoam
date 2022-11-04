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

#include "RK6TimeIntegrator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace timeIntegrators
{
    defineTypeNameAndDebug(RK6, 0);
    addToRunTimeSelectionTable(timeIntegrator, RK6, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrators::RK6::RK6
(
    const fvMesh& mesh,
    Istream& is
)
:
    timeIntegrator(mesh)
{
    // 3rd Choice: c2=c3=1/4, c4 = 2/4, c5=c6=3/4, c7=4/4
    this->as_ =
    {
        {1.0},
        {1.0, 0.0},
        {1.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
    };
    this->bs_ =
    {
        {0.25},
        {0.125, 0.125},
        {0.0, -5.0/6.0, 8.0/6.0},
        {0.125, 0.125, 0.0, 0.5},
        {0.0, 0.375, 0.25, -0.125, 0.25},
        {1.0/7.0, -2.0/7.0, 4.0/7.0, 0.0, 0.0, 4.0/7.0},
        {7.0/90.0, 0.0, 32.0/90.0, 12.0/90.0, 16.0/90.0, 16.0/90.0, 7.0/90.0}
    };
    this->initialize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrators::RK6::~RK6()
{}
// ************************************************************************* //
