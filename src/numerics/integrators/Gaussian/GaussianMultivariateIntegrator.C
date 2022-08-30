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

#include "GaussianIntegrator.H"
#include "GaussianWeights.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::GaussianMultivariateIntegrator<Type>::GaussianMultivariateIntegrator
(
    const equationType& eqn,
    const dictionary& dict
)
:
    MultivariateIntegrator<Type>(eqn, dict),
    ws_(eqn.nVar()),
    xs_(eqn.nVar())
{
    labelList nNodes(dict.lookup("nNodes"));
    forAll(nNodes, i)
    {
        GaussianWeights::setQuadrature(nNodes[i], ws_[i], xs_[i]);
    }
}


template<class Type>
Foam::GaussianMultivariateIntegrator<Type>::GaussianMultivariateIntegrator
(
    const equationType& eqn,
    const multivariateIntegrator& inter
)
:
    MultivariateIntegrator<Type>(eqn, inter),
    ws_(dynamicCast<const GaussianMultivariateIntegrator<Type>>(inter).ws_),
    xs_(dynamicCast<const GaussianMultivariateIntegrator<Type>>(inter).xs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::GaussianMultivariateIntegrator<Type>::integrateFunc
(
    const scalarList& x0,
    const scalarList& x1,
    const label li
) const
{
    labelList gi(x0.size(), 0);
    this->evals_++;
    Type fx(0.0*this->eqnPtr_->fX(x0, li));
    addCorners(0, x0, x1, li, gi, fx);
    return fx;
}


template<class Type>
void Foam::GaussianMultivariateIntegrator<Type>::addCorners
(
    const label diri,
    const scalarList& x0,
    const scalarList& x1,
    const label li,
    labelList& gi,
    Type& fx
) const
{
    if (diri < x0.size())
    {
        forAll(ws_[diri], i)
        {
            gi[diri] = i;
            addCorners(diri+1, x0, x1, li, gi, fx);
        }
        gi[diri] = 0;
        return;
    }

    scalarList x(x0);
    scalar w = 1.0;
    forAll(x0, i)
    {
        w *= ws_[i][gi[i]];
        x[i] = x0[i] + (x1[i] - x0[i])*xs_[i][gi[i]];
    }
    this->evals_++;
    fx = fx + w*this->eqnPtr_->fX(x, li);
}

template<class Type>
Type Foam::GaussianMultivariateIntegrator<Type>::integrate
(
    const scalarList& x0,
    const scalarList& x1,
    const label li
) const
{
    scalarList dX(x1 - x0);
    this->reset(dX);

    scalar dx = 1.0;
    forAll(x1, i)
    {
        dx *= x1[i] - x0[i];
    }
    Type Q(dx*integrateFunc(x0, x1, li));

    if (this->adaptive())
    {
        return this->integrate_(Q, 0, x0, x1, this->tolerance_, li);
    }
    return Q;
}

// ************************************************************************* //
