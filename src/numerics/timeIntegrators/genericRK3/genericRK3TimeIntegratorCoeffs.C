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

#include "genericRK3TimeIntegratorCoeffs.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace timeIntegrators
{
    defineTypeNameAndDebug(genericRK3, 0);
    addToRunTimeSelectionTable(timeIntegratorCoeffs, genericRK3, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrators::genericRK3::genericRK3(Istream& is)
:
    genericRK3(readScalar(is))
{}


Foam::timeIntegrators::genericRK3::genericRK3(const scalar alpha)
:
    timeIntegratorCoeffs(3),
    alpha_(alpha)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrators::genericRK3::~genericRK3()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeIntegrators::genericRK3::set
(
    List<List<scalar>>& as,
    List<List<scalar>>& bs,
    const label index
) const
{
    as = {{1.0}, {1.0, 0.0}, {1.0, 0.0, 0.0}};
    bs =
    {
        {alpha_},
        {
            1.0 + (1.0 - alpha_)/(alpha_*(3.0*alpha_ - 2.0)),
          - (1.0 - alpha_)/(alpha_*(3.0*alpha_ - 2.0))
        },
        {
            1.0 + 1.0/(6.0*alpha_),
            1.0/(6.0*alpha_*(1.0 - alpha_)),
            (2.0 - 3.0*alpha_)/((6.0*(1.0 - alpha_)))
        }
    };
}


// ************************************************************************* //
