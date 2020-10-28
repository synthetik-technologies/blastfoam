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

#include "RK3SSPTimeIntegrator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace timeIntegrators
{
    defineTypeNameAndDebug(RK3SSP, 0);
    addToRunTimeSelectionTable(timeIntegrator, RK3SSP, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrators::RK3SSP::RK3SSP
(
    const fvMesh& mesh,
    const label nSteps
)
:
    timeIntegrator(mesh, nSteps)
{
    if (nSteps <= 3)
    {
        this->as_ = {{1.0}, {0.75, 0.25}, {1.0/3.0, 0.0, 2.0/3.0}};
        this->bs_ = {{1.0}, {0.0, 0.25}, {0.0, 0.0, 2.0/3.0}};
    }
    else
    {
        if (nSteps > 4)
        {
            WarningInFunction
                << "RK3SSP only supports a maximum of 4 steps."
                << endl;
        }
        this->as_ =
        {
            {1.0},
            {0.0, 1.0},
            {2.0/3.0, 0.0, 1.0/3.0},
            {0.0, 0.0, 0.0, 1.0}
        };
        this->bs_ =
        {
            {0.5},
            {0.0, 0.5},
            {0.0, 0.0, 1.0/6.0},
            {0.0, 0.0, 0.0, 0.5}
        };
    }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrators::RK3SSP::~RK3SSP()
{}

// ************************************************************************* //
