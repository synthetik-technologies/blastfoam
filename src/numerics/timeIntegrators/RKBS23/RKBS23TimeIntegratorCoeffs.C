/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2023
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

#include "RKBS23TimeIntegratorCoeffs.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace timeIntegrators
{
    defineTypeNameAndDebug(RKBS23, 0);
    addToRunTimeSelectionTable(timeIntegratorCoeffs, RKBS23, dictionary);
    addToRunTimeSelectionTable
    (
        timeIntegratorCoeffs,
        RKBS23,
        dictionaryEmbedded
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrators::RKBS23::RKBS23(Istream& is)
:
    timeIntegratorCoeffs(4)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrators::RKBS23::~RKBS23()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeIntegrators::RKBS23::set
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
        {1.0, 0.0, 0.0, 0.0}
    };
    bs =
    {
        {1.0/2.0},
        {0.0, 3.0/4.0},
        {2.0/9.0, 1.0/3.0, 4.0/9.0},
        {2.0/9.0, 1.0/3.0, 4.0/9.0, 0.0},
        {7.0/24.0, 1.0/4.0, 1.0/3.0, 1.0/8.0}
    };
}
// ************************************************************************* //
