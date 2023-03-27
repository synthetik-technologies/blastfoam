/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2023
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

#include "RKF45TimeIntegratorCoeffs.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace timeIntegrators
{
    defineTypeNameAndDebug(RKF45, 0);
    addToRunTimeSelectionTable(timeIntegratorCoeffs, RKF45, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrators::RKF45::RKF45(Istream& is)
:
    timeIntegratorCoeffs(6)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrators::RKF45::~RKF45()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeIntegrators::RKF45::set
(
    List<List<scalar>>& as,
    List<List<scalar>>& bs,
    const label index
) const
{
    //- Fourth order coefficients
//     scalar b60 = 25.0/216.0;
//     scalar b61 = 0.0;
//     scalar b62 = 1408.0/2565.0;
//     scalar b63 = 2197.0/4104.0;
//     scalar b64 = -1.0/5.0;
//     scalar b65 = 0.0;

    //- Fifth order coefficients
    scalar b60 = 16.0/135.0;
    scalar b61 = 0.0;
    scalar b62 = 6656.0/12825.0;
    scalar b63 = 28561.0/56430.0;
    scalar b64 = -9.0/50.0;
    scalar b65 = 2.0/55.0;

    as =
    {
        {1.0},
        {1.0, 0.0},
        {1.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0}
    };
    bs =
    {
        {0.25},
        {3.0/32.0, 9.0/32.0},
        {1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0},
        {439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0},
        {-8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0},
        {b60, b61, b62, b63, b64, b65}
    };
}


// ************************************************************************* //
