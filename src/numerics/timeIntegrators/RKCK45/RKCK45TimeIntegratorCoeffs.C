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

#include "RKCK45TimeIntegratorCoeffs.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace timeIntegrators
{
    defineTypeNameAndDebug(RKCK45, 0);
    addToRunTimeSelectionTable(timeIntegratorCoeffs, RKCK45, dictionary);
    addToRunTimeSelectionTable
    (
        timeIntegratorCoeffs,
        RKCK45,
        dictionaryEmbedded
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrators::RKCK45::RKCK45(Istream& is)
:
    timeIntegratorCoeffs(6)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrators::RKCK45::~RKCK45()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeIntegrators::RKCK45::set
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
        {1.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0}
    };
    bs =
    {
        {1.0/5.0},
        {3.0/40.0, 9.0/40.0},
        {3.0/10.0, -9.0/10.0, 6.0/5.0},
        {-11.0/54.0, 5.0/2.0, -70.0/27.0, 35.0/27.0},
        {1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0},
        {37.0/378.0, 0.0, 250.0/621.0, 125.0/594.0, 0.0, 512.0/1771.0},
        {2825.0/27648.0, 0.0, 18575.0/48384.0, 13525.0/55296.0, 277.0/14336.0, 1.0/4.0}
    };
}


// ************************************************************************* //
