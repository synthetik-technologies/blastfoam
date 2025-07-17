/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "basisMinimizationScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace minimizationSchemes
{
    defineTypeNameAndDebug(basis, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::minimizationSchemes::basis::basis
(
    const scalarUnivariateEquation& eqns,
    const dictionary& dict
)
:
    minimizationScheme(eqns, dict),
    cmptLsEqn_(eqns, dict.subOrEmptyDict("basisLineSearchCoeffs"))
{}


Foam::minimizationSchemes::basis::basis
(
    const scalarUnivariateEquation& eqns,
    const basis& solver
)
:
    minimizationScheme(eqns, solver),
    cmptLsEqn_(eqns, solver.cmptLsEqn_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::minimizationSchemes::basis::searchDir
(
    const scalarList& x0,
    const label li,
    const label dir,
    const label sign,
    scalarList& xNew
) const
{
    cmptLsEqn_.setX0(x0);
    cmptLsEqn_.setCmpt(dir, sign);
    cmptLsEqn_.search(li, xNew);
}

// ************************************************************************* //
