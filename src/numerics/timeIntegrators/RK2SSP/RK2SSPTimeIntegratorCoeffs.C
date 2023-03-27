/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2023
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

#include "RK2SSPTimeIntegratorCoeffs.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace timeIntegrators
{
    defineTypeNameAndDebug(RK2SSP, 0);
    addToRunTimeSelectionTable(timeIntegratorCoeffs, RK2SSP, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrators::RK2SSP::RK2SSP(Istream& is)
:
    timeIntegratorCoeffs
    (
        !is.eof()
      ? readLabel(is)
      : 2
    )
{
    if (nSteps_ < 2)
    {
        WarningInFunction
            << "RK2SSP only supports a minimum of 2 steps."
            << endl;
        nSteps_ = 2;
    }
    else if (nSteps_ > 4)
    {
        WarningInFunction
            << "RK2SSP only supports a maximum of 4 steps."
            << endl;
        nSteps_ = 4;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrators::RK2SSP::~RK2SSP()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeIntegrators::RK2SSP::set
(
    List<List<scalar>>& as,
    List<List<scalar>>& bs,
    const label index
) const
{
    if (nSteps_ == 2)
    {
        as = {{1.0}, {0.5, 0.5}};
        bs = {{1.0}, {0.0, 0.5}};
    }
    else if (nSteps_ == 3)
    {
        as = {{1.0}, {0.0, 1.0}, {1.0/3.0, 0.0, 2.0/3.0}};
        bs = {{0.5}, {0.0, 0.5}, {0.0, 0.0, 1.0/3.0}};
    }
    else
    {
        as =
        {
            {1.0},
            {0.0, 1.0},
            {0.0, 0.0, 1.0},
            {0.25, 0.0, 0.0, 0.75}
        };
        bs =
        {
            {1.0/3.0},
            {0.0, 1.0/3.0},
            {0.0, 0.0, 1.0/3.0},
            {0.0, 0.0, 0.0, 0.25}
        };
    }
}


// ************************************************************************* //
