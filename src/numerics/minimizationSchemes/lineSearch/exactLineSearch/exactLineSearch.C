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

#include "exactLineSearch.H"
#include "univariateMinimizationScheme.H"
#include "goldenRatioUnivariateMinimizationScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(exactLineSearch, 0);
    addToRunTimeSelectionTable(lineSearch, exactLineSearch, dictionary);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::exactLineSearch::bracket
(
    const label li,
    const scalar x0
) const
{
    bracket(lsEqn_, x0, li, dx0_, k_, maxBracketIter_);
}


void Foam::exactLineSearch::bracket
(
    lineSearchEquation& eqn,
    const scalar x0,
    const label li,
    const scalar dx0,
    const scalar k,
    const label maxIter
)
{
    scalar dx = dx0;

    scalar xa = x0;
    scalar fxa = eqn.fx(xa, li);

    scalar xb = xa + dx;
    scalar fxb = eqn.fx(xb, li);

    if (mag(fxb - fxa) < small)
    {
        return;
    }

    if (fxb > fxa)
    {
        Swap(xa, xb);
        Swap(fxa, fxb);
        dx = -dx;
    }

    scalar xc = xb;
    scalar fxc = fxb;

    label iter = 0;
    while (iter++ < maxIter)
    {
        if (lineSearch::debug)
        {
            Info<< "Line search iteration " << iter << " "
                << ", direction = " << eqn.dir() << endl;
        }
        xc = xb + dx;
        fxc = eqn.fx(xc, li);

        if (fxc > fxb)
        {
            // Bound equation and return
            eqn.setLower(max(min(xa, xc), 0.0));
            eqn.setUpper(max(xa, xc));
            return;
        }

        xa = xb;
        fxa = fxb;

        xb = xc;
        fxb = fxc;

        dx *= k;
    }
    WarningInFunction
        << "Could not determine a valid interval containing a minimum in "
        << iter << " iterations " << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::exactLineSearch::exactLineSearch
(
    const scalarUnivariateEquation& eqns,
    const dictionary& dict
)
:
    lineSearch(eqns, dict),
    lineSearcher_
    (
        univariateMinimizationScheme::New
        (
            dict.lookupOrDefault
            (
                "solver",
                univariateMinimizationSchemes::goldenRatio::typeName
            ),
            lsEqn_,
            dict
        )
    ),
    maxBracketIter_(dict.lookupOrDefault("maxBracketIter", 100)),
    dx0_(dict.lookupOrDefault("dx0Bracket", 1e-2)),
    k_(dict.lookupOrDefault("kBracket", 2.0))
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::exactLineSearch::~exactLineSearch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::exactLineSearch::search
(
    const scalarList& x0,
    const scalarList& grad,
    const label li,
    scalarList& xNew
) const
{
    DebugInfo<< "Conducting line search" << endl;

    // Store current state of the line search class since this would
    // print out alot of information
    // Only print if debug level is sufficiently high
    const label oldDebug = univariateMinimizationScheme::debug;
    univariateMinimizationScheme::debug = lineSearch::debug;

    lsEqn_.update(x0, grad);
    bracket(li, 0.0);
    scalar alpha = lineSearcher_->solve(li);
    xNew = lsEqn_.calcX(alpha);

    // Reset debug flag
    univariateMinimizationScheme::debug = oldDebug;

    return;
}


void Foam::exactLineSearch::searchDir
(
    const scalarList& x0,
    const scalarList& dir,
    const label li,
    scalarList& xNew
) const
{
    DebugInfo<< "Conducting line search" << endl;

    // Store current state of the line search class since this would
    // print out alot of information
    // Only print if debug level is sufficiently high
    const label oldDebug = univariateMinimizationScheme::debug;
    univariateMinimizationScheme::debug = lineSearch::debug;

    lsEqn_.updateDir(x0, dir);
    scalar alpha = lineSearcher_->solve(li);
    xNew = lsEqn_.calcX(alpha);

    // Reset debug flag
    univariateMinimizationScheme::debug = oldDebug;
}


// ************************************************************************* //
