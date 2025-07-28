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

#include "RK4SSPTimeIntegratorCoeffs.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace timeIntegrators
{
    defineTypeNameAndDebug(RK4SSP, 0);
    addToRunTimeSelectionTable(timeIntegratorCoeffs, RK4SSP, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrators::RK4SSP::RK4SSP(Istream& is)
:
    timeIntegratorCoeffs(4)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrators::RK4SSP::~RK4SSP()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeIntegrators::RK4SSP::set
(
    List<List<scalar>>& as,
    List<List<scalar>>& bs,
    const label index
) const
{
    as =
    {
        {1.0},
        {649.0/1600.0, 951.0/1600.0},
        {53989.0/2500000.0, 4806213.0/20000000.0, 23619.0/32000.0},
        {0.2, 6127.0/30000.0, 7873.0/30000.0, 1.0/3.0}
    };
    bs =
    {
        {0.5},
        {-10890423.0/25193600.0, 5000.0/7873.0},
        {-102261.0/5000000.0, -5121.0/20000.0, 7873.0/10000.0},
        {0.1, 1.0/6.0, 0.0, 1.0/6.0}
    };
}


// ************************************************************************* //
