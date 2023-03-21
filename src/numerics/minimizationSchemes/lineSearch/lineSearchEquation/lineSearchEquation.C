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

#include "lineSearchEquation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char*
Foam::NamedEnum<Foam::lineSearchEquation::DescentMethod, 3>::names[] =
{
    "none",
    "FletcherReeves",
    "PolakRibiere"
};

const Foam::NamedEnum<Foam::lineSearchEquation::DescentMethod, 3>
Foam::lineSearchEquation::DescentMethodNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lineSearchEquation::lineSearchEquation
(
    const scalarUnivariateEquation& eqns,
    const word& descentType
)
:
    ScalarEquation(-great, great),
    eqns_(eqns),
    descent_(DescentMethodNames_[descentType]),
    grad_(eqns_.nVar(), 0.0),
    dir_(eqns_.nVar(), 0.0),
    x0_(eqns_.lowerLimits()),
    beta_(1.0),
    x_(x0_),
    dfdX_(eqns_.nVar(), 0.0)
{
    if (isA<scalarEquation>(eqns))
    {
        sEqn_.set
        (
            &dynamicCast<const scalarEquation>(eqns)
        );
    }
}

Foam::lineSearchEquation::lineSearchEquation
(
    const scalarUnivariateEquation& eqns,
    const DescentMethod descent
)
:
    ScalarEquation(-great, great),
    eqns_(eqns),
    descent_(descent),
    grad_(eqns_.nVar(), 0.0),
    dir_(eqns_.nVar(), 0.0),
    x0_(eqns_.lowerLimits()),
    beta_(1.0),
    x_(x0_),
    dfdX_(eqns_.nVar(), 0.0)
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

Foam::lineSearchEquation::~lineSearchEquation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lineSearchEquation::update
(
    const scalarList& x0,
    const scalarList& grad
)
{
    // Set reference point
    x0_ = x0;

    // Set travel direction, aka normalized, negative gradient
    if (descent_ == NONE)
    {
        scalar gradMagSqr = 0.0;
        forAll(grad_, i)
        {
            gradMagSqr += sqr(grad[i]);
        }
        dir_ = grad;
        dir_ /= -max(sqrt(gradMagSqr), small);
        beta_ = 1.0;
    }
    else
    {
        scalar gradOldMagSqr = 0.0;
        forAll(grad_, i)
        {
            gradOldMagSqr += sqr(grad_[i]);
        }
        if (descent_ == FLETCHER_REEVES)
        {
            scalar gradMagSqr = 0.0;
            forAll(grad, i)
            {
                gradMagSqr += sqr(grad[i]);
            }
            beta_ = gradMagSqr/max(gradOldMagSqr, small);
        }
        else if (descent_ == POLAK_RIBIERE)
        {
            beta_ = 0.0;
            forAll(grad_, i)
            {
                beta_ += grad[i]*(grad[i] - grad_[i]);
            }
            beta_ /= max(gradOldMagSqr, small);
        }

        // Update the travel direction
        forAll(dir_, i)
        {
            dir_[i] = -grad[i] + beta_*dir_[i];
        }
    }

    // Compute total distance, and distances to the true bounds
    scalar maxDistSqr = 0.0;
    forAll(grad, i)
    {
        maxDistSqr +=
            max
            (
                sqr((eqns_.upperLimits()()[i] - x0[i])*dir_[i]),
                sqr((eqns_.lowerLimits()()[i] - x0[i])*dir_[i])
            );
    }

    // Normal component of the distance to the boundary
    this->setUpper(sqrt(maxDistSqr));

    // Store the gradient
    grad_ = grad;
}


void Foam::lineSearchEquation::updateDir
(
    const scalarList& x0,
    const scalarList& dir
)
{
    x0_ = x0;
    dir_ = dir;

    // Compute total distance, and distances to the true bounds
    tmp<scalarField> lowerLimits(eqns_.lowerLimits());
    tmp<scalarField> upperLimits(eqns_.upperLimits());
    scalar lower = 0.0;
    scalar upper = 0.0;
    forAll(dir_, i)
    {
        lower += (lowerLimits()[i] - x0[i])*dir_[i];
        upper += (upperLimits()[i] - x0[i])*dir_[i];
    }

    // if (lower > upper)
    // {
    //     Swap(lower, upper);
    // }

    // Normal component of the distance to the boundary
    this->setLower(lower);
    this->setUpper(upper);
}


const Foam::scalarField& Foam::lineSearchEquation::calcX
(
    const scalar dist
) const
{
    forAll(x_, i)
    {
        x_[i] = x0_[i] + dist*dir_[i];
    }
    return x_;
}


// ************************************************************************* //
