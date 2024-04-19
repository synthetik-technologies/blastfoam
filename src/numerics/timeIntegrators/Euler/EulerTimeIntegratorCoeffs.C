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

#include "EulerTimeIntegratorCoeffs.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace timeIntegrators
{
    defineTypeNameAndDebug(Euler, 0);
    addToRunTimeSelectionTable(timeIntegratorCoeffs, Euler, dictionary);
    addToRunTimeSelectionTable
    (
        timeIntegratorCoeffs,
        Euler,
        dictionaryEmbedded
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrators::Euler::Euler(Istream& is)
:
    timeIntegratorCoeffs(1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrators::Euler::~Euler()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeIntegrators::Euler::set
(
    List<List<scalar>>& as,
    List<List<scalar>>& bs,
    const label index
) const
{
    as = {{1.0}};
    bs = {{1.0}, {0.0}};
}

// ************************************************************************* //
