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

#include "AdamsBashforth4TimeIntegratorCoeffs.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace timeIntegrators
{
    defineTypeNameAndDebug(AdamsBashforth4, 0);
    addToRunTimeSelectionTable
    (
        timeIntegratorCoeffs,
        AdamsBashforth4,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrators::AdamsBashforth4::AdamsBashforth4(Istream& is)
:
    AdamsBashforthBase(4)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrators::AdamsBashforth4::~AdamsBashforth4()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeIntegrators::AdamsBashforth4::set
(
    List<List<scalar>>& as,
    List<List<scalar>>& bs,
    const label index
) const
{
    if (index <= 0)
    {
        as = {{1.0}};
        bs = {{1.0}};
    }
    else if (index == 1)
    {
        as = {{1.0, 0.0}};
        bs = {{1.5, -0.5}};
    }
    else if (index == 2)
    {
        scalar s = 1.0/12.0;
        as = {{1.0, 0.0, 0.0}};
        bs = {{23.0*s, -16.0*s, 5.0*s}};
    }
    else //if (index == 3)
    {
        scalar s = 1.0/24.0;
        as = {{1.0, 0.0, 0.0, 0.0}};
        bs = {{55.0*s, -59.0*s, 37.0*s, -9.0*s}};
    }
    // else
    // {
    //
    //     as = {{1.0, 0.0, 0.0, 0.0, 0.0}};
    //     bs = {{251.0*s, 646.0*s, -264.0*s, 106.0*s, -19.0*s}};
    // }
}


// ************************************************************************* //
