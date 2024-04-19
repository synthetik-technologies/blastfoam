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

#include "RKF12TimeIntegratorCoeffs.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace timeIntegrators
{
    defineTypeNameAndDebug(RKF12, 0);
    addToRunTimeSelectionTable(timeIntegratorCoeffs, RKF12, dictionary);
    addToRunTimeSelectionTable
    (
        timeIntegratorCoeffs,
        RKF12,
        dictionaryEmbedded
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrators::RKF12::RKF12(Istream& is)
:
    timeIntegratorCoeffs(3)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrators::RKF12::~RKF12()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeIntegrators::RKF12::set
(
    List<List<scalar>>& as,
    List<List<scalar>>& bs,
    const label index
) const
{
    as =
    {
        {1.0},
        {1.0, 0.0},
        {1.0, 0.0, 0.0},
    };
    bs =
    {
        {1.0/2.0},
        {1.0/256.0, 255.0/256.0},
        {1.0/512.0, 255/256.0, 1.0/512.0},
        {1.0/256.0, 255.0/256.0, 0.0}

    };
}


// ************************************************************************* //
