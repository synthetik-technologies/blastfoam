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

#include "GaussianIntegrator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
void Foam::GaussianIntegrator<Type>::setQuadrature(const label N)
{
    switch (N)
    {
        case 1:
        {
            xs_ = {0.0};
            ws_ = {2.0};
            break;
        }
        case 2:
        {
            const scalar x = 1.0/sqrt(3.0);
            xs_ = {-x, x};
            ws_ = {1.0, 1.0};
            break;
        }
        case 3:
        {
            const scalar x = sqrt(3.0/5.0);
            xs_ = {-x, 0.0,  x};
            ws_ = {5.0/9.0, 8.0/9.0, 5.0/9.0};
            break;
        }
        case 4:
        {
            const scalar x1 = sqrt(3.0/7.0 - 2.0/7.0*sqrt(6.0/5.0));
            const scalar x2 = sqrt(3.0/7.0 + 2.0/7.0*sqrt(6.0/5.0));

            const scalar w1 = (18.0 + sqrt(30.0))/36.0;
            const scalar w2 = (18.0 - sqrt(30.0))/36.0;

            xs_ = {-x2, -x1, x1, x2};
            ws_ = {w2, w1, w1, w2};
            break;
        }
        case 5:
        {
            const scalar x1 = sqrt(5.0 - 2.0*sqrt(10.0/7.0))/3.0;
            const scalar x2 = sqrt(5.0 + 2.0*sqrt(10.0/7.0))/3.0;

            const scalar w0 = 128.0/225.0;
            const scalar w1 = (322.0 + 13.0*sqrt(70.0))/900.0;
            const scalar w2 = (322.0 - 13.0*sqrt(70.0))/900.0;

            xs_ = {-x2, -x1, 0.0, x1, x2};
            ws_ = {w2, w1, w0, w1, w2};
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unsupported number of Nodes for Gaussian integration" << nl
                << "1 to 5 nodes can be used."
                << abort(FatalError);
        }
    }

    // Normalize weights and move xs to be in (0, 1)
    scalar sumW = sum(ws_);
    forAll(xs_, i)
    {
        ws_[i] /= sumW;
        xs_[i] = 0.5*(xs_[i] + 1.0);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::GaussianIntegrator<Type>::GaussianIntegrator
(
    const equationType& eqn,
    const dictionary& dict
)
:
    Integrator<Type>(eqn, dict),
    ws_(0),
    xs_(0)
{
    setQuadrature(dict.lookup<label>("nNodes"));
}


template<class Type>
Foam::GaussianIntegrator<Type>::GaussianIntegrator
(
    const equationType& eqn,
    const integrator& inter
)
:
    Integrator<Type>(eqn, inter),
    ws_(dynamicCast<const GaussianIntegrator<Type>>(inter).ws_),
    xs_(dynamicCast<const GaussianIntegrator<Type>>(inter).xs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::GaussianIntegrator<Type>::integrate_
(
    const scalar dx,
    const PtrList<Type>& fs
) const
{
    Type res(ws_[0]*fs[0]);
    for (label i = 1; i < ws_.size(); i++)
    {
        res = res + ws_[i]*fs[i];
    }
    return dx*res;
}


template<class Type>
Type Foam::GaussianIntegrator<Type>::integrate_
(
    const Type& Q,
    const scalar x0,
    const scalar x1,
    const scalarList& xs,
    const PtrList<Type>& fs,
    const scalar tol,
    const label li
) const
{
    const scalar dx(x1 - x0);
    if (mag(dx) < this->minDx_)
    {
        return Q;
    }

    const scalar xm = x0 + 0.5*dx;
    scalarList x0s(xs.size());
    scalarList x1s(xs.size());
    PtrList<Type> f0s(fs.size());
    PtrList<Type> f1s(fs.size());
    const scalar dxBy2 = dx/2.0;
    forAll(xs, i)
    {
        x0s[i] = x0 + dxBy2*xs_[i];
        x1s[i] = xm + dxBy2*xs_[i];

        f0s.set(i, new Type(this->eqnPtr_->fx(x0s[i], li)));
        f1s.set(i, new Type(this->eqnPtr_->fx(x1s[i], li)));
    }
    this->evals_ += 2.0*xs_.size();

    const Type fx0(integrate_(dxBy2, f0s));
    const Type fx1(integrate_(dxBy2, f1s));
    const Type fx(fx0 + fx1);

    if (this->converged(fx, Q, dx, tol))
    {
        return fx;
    }
    else
    {
        this->intervals_++;
        return
            integrate_(fx0, x0, xm, x0s, f0s, tol/2.0, li)
          + integrate_(fx1, xm, x1, x1s, f1s, tol/2.0, li);
    }
}


template<class Type>
Type Foam::GaussianIntegrator<Type>::integrate
(
    const scalar x0,
    const scalar x1,
    const label li
) const
{
    scalar dx(x1 - x0);
    if (mag(dx) < small)
    {
        return dx*this->eqnPtr_->fx(x0, li);
    }

    this->reset(dx);

    if (this->adaptive())
    {
        scalarList xs(xs_.size());
        PtrList<Type> fs(xs_.size());
        forAll(xs, i)
        {
            xs[i] = x0 + dx*xs_[i];
            fs.set(i, new Type(this->eqnPtr_->fx(xs[i], li)));
        }
        this->evals_ = xs_.size();

        return integrate_
        (
            integrate_(dx, fs),
            x0, x1,
            xs,
            fs,
            this->tolerance_,
            li
        );
    }

    dx /= scalar(this->nIntervals_);
    scalarList xs(xs_.size());
    PtrList<Type> fs(xs_.size());
    forAll(xs, i)
    {
        xs[i] = x0 + dx*xs_[i];
        fs.set(i, new Type(this->eqnPtr_->fx(xs[i], li)));
    }

    Type res(integrate_(dx, fs));
    for (label i = 1; i < this->nIntervals_; i++)
    {
        forAll(xs, i)
        {
            xs[i] += dx;
            fs[i] = this->eqnPtr_->fx(xs[i], li);
        }

        res = res + integrate_(dx, fs);
    }
    this->evals_ = xs_.size()*this->nIntervals_;
    return res;
}

// ************************************************************************* //
