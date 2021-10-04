/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
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

#include "NewtonRaphsonMinimizationScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(NewtonRaphsonMinimizationScheme, 0);
    addToRunTimeSelectionTable
    (
        minimizationScheme,
        NewtonRaphsonMinimizationScheme,
        dictionaryTwo
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NewtonRaphsonMinimizationScheme::NewtonRaphsonMinimizationScheme
(
    const scalarEquation& eqn,
    const dictionary& dict
)
:
    minimizationScheme(eqn, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::NewtonRaphsonMinimizationScheme::solve
(
    const scalar x,
    const scalar xLow,
    const scalar xHigh,
    const label li
) const
{
    scalar xOld = x;
    scalar xNew = x;
    for (stepi_ = 0; stepi_ < maxSteps_; stepi_++)
    {
        xNew = xOld - eqn_.dfdx(xOld, li)/stabilise(eqn_.d2fdx2(xOld, li), small);

        error_ = mag(xNew - xOld);
        limit(xNew);
        if (error_ < tolerance_)
        {
            return xNew;
        }
        xOld = xNew;

    }
    WarningInFunction
        << "Could not converge to the given root." << endl;

    return xNew;
}

// ************************************************************************* //
