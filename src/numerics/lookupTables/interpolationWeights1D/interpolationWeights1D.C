/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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

#include "interpolationWeights1D.H"
#include "hashedWordList.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::interpolationWeight1D> Foam::interpolationWeight1D::New
(
    const word& scheme,
    const List<scalar>& xs
)
{
    if (scheme == "floor")
    {
        return autoPtr<interpolationWeight1D>
        (
            new interpolationWeights1D::floor(xs)
        );
    }
    else if (scheme == "ceil")
    {
        return autoPtr<interpolationWeight1D>
        (
            new interpolationWeights1D::ceil(xs)
        );
    }
    else if (scheme == "linearClamp")
    {
        return autoPtr<interpolationWeight1D>
        (
            new interpolationWeights1D::linearClamp(xs)
        );
    }
    else if (scheme == "linearExtrapolated")
    {
        return autoPtr<interpolationWeight1D>
        (
            new interpolationWeights1D::linearExtrapolated(xs)
        );
    }
    else if (scheme == "quadraticClamp")
    {
        return autoPtr<interpolationWeight1D>
        (
            new interpolationWeights1D::quadraticClamp(xs)
        );
    }
    else if (scheme == "quadraticExtrapolated")
    {
        return autoPtr<interpolationWeight1D>
        (
            new interpolationWeights1D::quadraticExtrapolated(xs)
        );
    }
    else if (scheme == "cubicClamp")
    {
        return autoPtr<interpolationWeight1D>
        (
            new interpolationWeights1D::cubicClamp(xs)
        );
    }
    else if (scheme == "cubicExtrapolated")
    {
        return autoPtr<interpolationWeight1D>
        (
            new interpolationWeights1D::cubicExtrapolated(xs)
        );
    }
    else
    {
        FatalErrorInFunction
            << scheme << " is not a valid interpolation scheme" << nl
            << "Options are: " << nl
            << "    linearClamp" << nl
            << "    linearExtrapolated" << nl
            << "    quadraticClamp" << nl
            << "    quadraticExtrapolated" << nl
            << "    cubicClamp" << nl
            << "    cubicExtrapolated" << nl
            << "    ceil" << nl
            << "    floor" << nl
            << abort(FatalError);
    }
    return autoPtr<interpolationWeight1D>();
}

bool Foam::interpolationWeight1D::validate(const bool fail) const
{
    hashedWordList validSchemes;
    const label n = xs_.size();
    if (n > 0)
    {
        validSchemes.append("floor");
        validSchemes.append("ceil");
        validSchemes.append("linearClamp");
    }
    if (n > 1)
    {
        validSchemes.append("linearExtrapolated");
    }
    if (n > 2)
    {
        validSchemes.append("quadraticClamp");
        validSchemes.append("quadraticExtrapolated");
    }
    if (n > 3)
    {
        validSchemes.append("cubicClamp");
        validSchemes.append("cibicExtrapolated");
    }

    if (!validSchemes.found(this->type()))
    {
        if (fail)
        {
            FatalErrorInFunction
                << this->type()
                << " is not a valid interpolation scheme for a list "
                << " of size " << n << nl
                << "Valid Options are: " << nl
                << validSchemes << endl
                << abort(FatalError);
        }
        return false;
    }
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::interpolationWeights1D::floor::updateWeights
(
    const scalar x,
    const label i,
    List<label>& indices,
    List<scalar>& weights
) const
{
    indices.setSize(1);
    weights.setSize(1);

    indices[0] = i;
    weights[0] = 1.0;
}


void Foam::interpolationWeights1D::ceil::updateWeights
(
    const scalar x,
    const label i,
    List<label>& indices,
    List<scalar>& weights
) const
{
    indices.setSize(1);
    weights.setSize(1);

    indices[0] = x != xs_[i] ? i + 1 : i;
    weights[0] = 1.0;
}


void Foam::interpolationWeights1D::linearExtrapolated::updateWeights
(
    const scalar x,
    const label i,
    List<label>& indices,
    List<scalar>& weights
) const
{
    label lo = max(i, 0);
    label hi = min(xs_.size() - 1, lo + 1);

    indices.setSize(2);
    weights.setSize(2);

    indices[0] = lo;
    indices[1] = hi;

    weights[1] = (x - xs_[lo])/(xs_[hi] - xs_[lo]);
    weights[0] = 1.0 - weights[1];
}


void Foam::interpolationWeights1D::quadraticExtrapolated::updateWeights
(
    const scalar x,
    const label i,
    List<label>& indices,
    List<scalar>& weights
) const
{
    label lo = max(i-1, 0);
    label mid = lo + 1;
    label hi = mid + 1;

    const scalar& x0 = xs_[lo];
    const scalar& x1 = xs_[mid];
    const scalar& x2 = xs_[hi];

    indices.resize(3);
    weights.resize(3);

    indices[0] = lo;
    indices[1] = mid;
    indices[2] = hi;

    weights[0] = (x - x1)*(x - x2)/(x0 - x1)/(x0 - x2);
    weights[1] = (x - x2)*(x - x0)/(x1 - x2)/(x1 - x0);
    weights[2] = (x - x0)*(x - x1)/(x2 - x0)/(x2 - x1);
}


void Foam::interpolationWeights1D::cubicExtrapolated::updateWeights
(
    const scalar x,
    const label i,
    List<label>& indices,
    List<scalar>& weights
) const
{
    label lo = max(i - 1, 0);
    label hi = lo + 3;
    if (hi > xs_.size() - 1)
    {
        hi = xs_.size() - 1;
        lo = hi - 3;
    }

    indices.resize(4);
    weights.resize(4);

    indices[0] = lo;
    indices[1] = lo + 1;
    indices[2] = lo + 2;
    indices[3] = lo + 3;

    const scalar& x0 = xs_[lo];
    const scalar& x1 = xs_[lo+1];
    const scalar& x2 = xs_[lo+2];
    const scalar& x3 = xs_[lo+3];


    weights[0] = (x - x1)*(x - x2)*(x - x3)/(x0 - x1)/(x0 - x2)/(x0 - x3);
    weights[1] = (x - x2)*(x - x3)*(x - x0)/(x1 - x0)/(x1 - x2)/(x1 - x3);
    weights[2] = (x - x3)*(x - x0)*(x - x1)/(x2 - x0)/(x2 - x1)/(x2 - x3);
    weights[3] = (x - x0)*(x - x1)*(x - x2)/(x3 - x0)/(x3 - x1)/(x3 - x2);

//     const scalar dx20 = xs[I+1] - xs[I-1];
//     const scalar dx31 = xs[I+2] - xs[I];
//
//     const scalar x3 = pow3(x);
//     const scalar x2 = sqr(x);
//
//     weights[0] = (-x3 + x2*2.0 - x)/dx20;
//     weights[1] = x3*(2.0 - 1.0/dx31) + x2*(1.0/dx31 - 3.0) + 1.0;
//     weights[2] = x3*(1.0/dx20 - 2.0) + x2*(3.0 - 2.0/dx20) + 1.0/dx20;
//     weights[3] = (x3 - x2)/dx31;
}

// ************************************************************************* //
