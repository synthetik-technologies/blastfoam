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

#include "RK1SSPTimeIntegratorCoeffs.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace timeIntegrators
{
    defineTypeNameAndDebug(RK1SSP, 0);
    addToRunTimeSelectionTable(timeIntegratorCoeffs, RK1SSP, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrators::RK1SSP::RK1SSP(Istream& is)
:
    timeIntegratorCoeffs
    (
        !is.eof()
      ? readLabel(is)
      : 1
    )
{
    if (nSteps_ < 1)
    {
        nSteps_ = 1;
    }
    if (nSteps_ > 3)
    {
        WarningInFunction
            << "RK1SSP only supports a maximum of 3 steps."
            << endl;
        nSteps_ = 3;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrators::RK1SSP::~RK1SSP()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeIntegrators::RK1SSP::set
(
    List<List<scalar>>& as,
    List<List<scalar>>& bs,
    const label index
) const
{
    if (nSteps_ == 1)
    {
        as = {{1.0}};
        bs = {{1.0}};
    }
    else if (nSteps_ == 2)
    {
        as = {{1.0}, {0.0, 1.0}};
        bs = {{0.5}, {0.0, 0.5}};
    }
    else
    {
        as =
        {
            {1.0},
            {0.0, 1.0},
            {0.0, 0.0, 1.0},
        };
        bs =
        {
            {1.0/3.0},
            {0.0, 1.0/3.0},
            {0.0, 0.0, 1.0/3.0},
        };
    }
}


// ************************************************************************* //
