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
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(interpolationWeight1D, 0);
    defineRunTimeSelectionTable(interpolationWeight1D, null);
    defineRunTimeSelectionTable(interpolationWeight1D, dictionary);
}

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::interpolationWeight1D> Foam::interpolationWeight1D::New
(
    const word& scheme,
    const List<scalar>& xs,
    const bool finalData
)
{
    DebugInfo << "Selecting interpolation scheme: " << scheme << endl;

    nullConstructorTable::iterator cstrIter =
        nullConstructorTablePtr_->find(scheme);

    if (cstrIter == nullConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown interpolation scheme " << scheme << nl << nl
            << "Valid interpolation schemes are: " << endl
            << nullConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return autoPtr<interpolationWeight1D>(cstrIter()(xs));
}


Foam::autoPtr<Foam::interpolationWeight1D> Foam::interpolationWeight1D::New
(
    const word& scheme,
    const dictionary& dict,
    const List<scalar>& xs,
    const bool finalData
)
{
    DebugInfo << "Selecting interpolation scheme: " << scheme << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(scheme);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown interpolation scheme " << scheme << nl << nl
            << "Valid interpolation schemes are: " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return autoPtr<interpolationWeight1D>(cstrIter()(dict, xs));
}


bool Foam::interpolationWeight1D::validate(const bool fail) const
{
    hashedWordList validSchemes;
    const label n = xs_.size();
    if (n >= nRequired())
    {
        return true;
    }
    else if (!fail)
    {
        return false;
    }

    forAllConstIter(nullConstructorTable, *nullConstructorTablePtr_, iter)
    {
        autoPtr<interpolationWeight1D> scheme(iter()(xs_));
        if (scheme->nRequired() <= n)
        {
            validSchemes.append(scheme->type());
        }
    }
    FatalErrorInFunction
        << this->type()
        << " requires " << this->nRequired() << " entries, but "
        << n << " entries were provied" << nl << nl
        << "Valid schemes for " << n << " entries are: " << nl
        << validSchemes << endl
        << abort(FatalError);
    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace interpolationWeights1D
{
    defineTypeNameAndDebug(floor, 0);
    addToRunTimeSelectionTable(interpolationWeight1D, floor, null);
    addToRunTimeSelectionTable(interpolationWeight1D, floor, dictionary);
}
}

void Foam::interpolationWeights1D::floor::updateWeights
(
    const scalar x,
    const label i,
    DynamicList<label>& indices,
    DynamicList<scalar>& weights
) const
{
    indices.setSize(1);
    weights.setSize(1);

    indices[0] = i;
    weights[0] = 1.0;
}


namespace Foam
{
namespace interpolationWeights1D
{
    defineTypeNameAndDebug(ceil, 0);
    addToRunTimeSelectionTable(interpolationWeight1D, ceil, null);
    addToRunTimeSelectionTable(interpolationWeight1D, ceil, dictionary);
}
}

void Foam::interpolationWeights1D::ceil::updateWeights
(
    const scalar x,
    const label i,
    DynamicList<label>& indices,
    DynamicList<scalar>& weights
) const
{
    indices.setSize(1);
    weights.setSize(1);

    indices[0] = x != xs_[i] ? i + 1 : i;
    weights[0] = 1.0;
}


namespace Foam
{
namespace interpolationWeights1D
{
    defineTypeNameAndDebug(linearExtrapolated, 0);
    defineTypeNameAndDebug(linearClamp, 0);
    addToRunTimeSelectionTable
    (
        interpolationWeight1D,
        linearExtrapolated,
        null
    );
    addToRunTimeSelectionTable
    (
        interpolationWeight1D,
        linearExtrapolated,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        interpolationWeight1D,
        linearClamp,
        null
    );
    addToRunTimeSelectionTable
    (
        interpolationWeight1D,
        linearClamp,
        dictionary
    );
}
}


void Foam::interpolationWeights1D::linearExtrapolated::updateWeights
(
    const scalar x,
    const label i,
    DynamicList<label>& indices,
    DynamicList<scalar>& weights
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


namespace Foam
{
namespace interpolationWeights1D
{
    defineTypeNameAndDebug(quadraticExtrapolated, 0);
    defineTypeNameAndDebug(quadraticClamp, 0);
    addToRunTimeSelectionTable
    (
        interpolationWeight1D,
        quadraticExtrapolated,
        null
    );
    addToRunTimeSelectionTable
    (
        interpolationWeight1D,
        quadraticExtrapolated,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        interpolationWeight1D,
        quadraticClamp,
        null
    );
    addToRunTimeSelectionTable
    (
        interpolationWeight1D,
        quadraticClamp,
        dictionary
    );
}
}


void Foam::interpolationWeights1D::quadraticExtrapolated::updateWeights
(
    const scalar x,
    const label i,
    DynamicList<label>& indices,
    DynamicList<scalar>& weights
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


namespace Foam
{
namespace interpolationWeights1D
{
    defineTypeNameAndDebug(cubicExtrapolated, 0);
    defineTypeNameAndDebug(cubicClamp, 0);
    addToRunTimeSelectionTable
    (
        interpolationWeight1D,
        cubicExtrapolated,
        null
    );
    addToRunTimeSelectionTable
    (
        interpolationWeight1D,
        cubicExtrapolated,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        interpolationWeight1D,
        cubicClamp,
        null
    );
    addToRunTimeSelectionTable
    (
        interpolationWeight1D,
        cubicClamp,
        dictionary
    );
}
}

void Foam::interpolationWeights1D::cubicExtrapolated::updateWeights
(
    const scalar x,
    const label i,
    DynamicList<label>& indices,
    DynamicList<scalar>& weights
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
