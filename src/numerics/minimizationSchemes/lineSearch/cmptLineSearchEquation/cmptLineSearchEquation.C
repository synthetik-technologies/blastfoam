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

#include "cmptLineSearchEquation.H"
#include "univariateMinimizationScheme.H"
#include "goldenRatioUnivariateMinimizationScheme.H"
#include "exactLineSearch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cmptLineSearchEquation::cmptLineSearchEquation
(
    const scalarUnivariateEquation& eqns,
    const dictionary& dict
)
:
    ScalarEquation(-great, great),
    eqns_(eqns),
    cmpt_(0),
    sign_(1),
    x0_(eqns_.lowerLimits()),
    x_(x0_),
    dfdX_(eqns_.nVar(), 0.0),
    lineSearcher_
    (
        univariateMinimizationScheme::New
        (
            dict.lookupOrDefault
            (
                "solver",
                univariateMinimizationSchemes::goldenRatio::typeName
            ),
            *this,
            dict
        )
    )
{
    if (isA<scalarEquation>(eqns))
    {
        sEqn_.set
        (
            &dynamicCast<const scalarEquation>(eqns)
        );
    }
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::cmptLineSearchEquation::~cmptLineSearchEquation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cmptLineSearchEquation::setX0(const scalarList& x0)
{
    x0_ = x0;
}


void Foam::cmptLineSearchEquation::setCmpt
(
    const label cmpt,
    const label sign
)
{
    cmpt_ = cmpt;
    sign_ = sign;

    this->setLower(eqns_.lowerLimits()()[cmpt_] - x0_[cmpt_]);
    this->setUpper(eqns_.upperLimits()()[cmpt_] - x0_[cmpt_]);
}


void Foam::cmptLineSearchEquation::flip()
{
    sign_ = -sign_;
    this->setLower(this->upper());
    this->setUpper(this->lower());
}


const Foam::scalarField& Foam::cmptLineSearchEquation::calcX
(
    const scalar dist
) const
{
    x_ = x0_;
    x_[cmpt_] = x0_[cmpt_] + dist*sign_;
    return x_;
}


void Foam::cmptLineSearchEquation::search
(
    const label li,
    scalarList& xNew
)
{
    DebugInfo<< "Conducting line search" << endl;

    // Store current state of the line search class since this would
    // print out alot of information
    // Only print if debug level is sufficiently high
    const label oldDebug = univariateMinimizationScheme::debug;
    univariateMinimizationScheme::debug = lineSearch::debug;

    xNew = calcX(lineSearcher_->solve(li));

    // Reset debug flag
    univariateMinimizationScheme::debug = oldDebug;
}

// ************************************************************************* //
