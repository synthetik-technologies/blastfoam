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

#include "genericRK3TimeIntegrator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace timeIntegrators
{
    defineTypeNameAndDebug(genericRK3, 0);
    addToRunTimeSelectionTable(timeIntegrator, genericRK3, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrators::genericRK3::genericRK3
(
    const fvMesh& mesh,
    Istream& is
)
:
    timeIntegrator(mesh)
{
    scalar alpha(readScalar(is));

    this->as_ = {{1.0}, {1.0, 0.0}, {1.0, 0.0, 0.0}};
    this->bs_ =
    {
        {alpha},
        {
            1.0 + (1.0 - alpha)/(alpha*(3.0*alpha - 2.0)),
          - (1.0 - alpha)/(alpha*(3.0*alpha - 2.0))
        },
        {
            1.0 + 1.0/(6.0*alpha),
            1.0/(6.0*alpha*(1.0 - alpha)),
            (2.0 - 3.0*alpha)/((6.0*(1.0 - alpha)))
        }
    };
    initialize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrators::genericRK3::~genericRK3()
{}
// ************************************************************************* //
