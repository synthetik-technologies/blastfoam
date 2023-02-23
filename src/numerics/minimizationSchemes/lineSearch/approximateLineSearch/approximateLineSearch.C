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

#include "approximateLineSearch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(approximateLineSearch, 0);
    addToRunTimeSelectionTable(lineSearch, approximateLineSearch, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::approximateLineSearch::approximateLineSearch
(
    const scalarUnivariateEquation& eqns,
    const dictionary& dict
)
:
    lineSearch(eqns, dict),
    alpha0_(dict.lookupOrDefault("alpha0", 1.0)),
    alpha_(alpha0_),
    beta_(dict.lookupOrDefault("beta", 1e-4)),
    p_(dict.lookupOrDefault("p", 0.5))
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::approximateLineSearch::~approximateLineSearch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::approximateLineSearch::search
(
    const scalarList& x0,
    const scalarList& grad,
    const label li,
    scalarList& xNew
) const
{
    if (debug)
    {
        Info<< "Conducting approximate line search" << endl;
    }

    scalar fx0 = eqns_.fX(x0, li);
    alpha_ = alpha0_;

    lsEqn_.update(x0, grad);

    scalar dfx = 0.0;
    forAll(grad, i)
    {
        dfx += lsEqn_.dir()[i]*grad[i];
    }

    while (lsEqn_.fx(alpha_, li) > fx0 + alpha_*beta_*dfx)
    {
        alpha_ *= p_;
    }

    if (alpha_ < small)
    {
        alpha_ = -alpha0_;
        while (lsEqn_.fx(alpha_, li) > fx0 + alpha_*beta_*dfx)
        {
            alpha_ *= p_;
        }
    }

    Info<<"alpha: "<<alpha_<<" "<<grad<<" "<<lsEqn_.dir()<<endl;
    xNew = lsEqn_.calcX(alpha_);
}


// ************************************************************************* //
