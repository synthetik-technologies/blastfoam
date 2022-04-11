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

#include "MultivariateIntegrator.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::MultivariateIntegrator<Type>>
Foam::MultivariateIntegrator<Type>::New
(
    const equationType& eqn,
    const dictionary& dict
)
{
    return New(dict.lookup<word>("integrator"), eqn, dict);
}


template<class Type>
Foam::autoPtr<Foam::MultivariateIntegrator<Type>>
Foam::MultivariateIntegrator<Type>::New
(
    const word& integratorTypeName,
    const equationType& eqn,
    const dictionary& dict
)
{
    Info<< "Selecting integrator " << integratorTypeName << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(integratorTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown integrator type "
            << integratorTypeName << nl << nl
            << "Valid integrators are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<MultivariateIntegrator<Type>>(cstrIter()(eqn, dict));
}


template<class Type>
Foam::autoPtr<Foam::MultivariateIntegrator<Type>> Foam::MultivariateIntegrator<Type>::New
(
    const equationType& eqn,
    const multivariateIntegrator& inter
)
{
    typename inputsConstructorTable::iterator cstrIter =
        inputsConstructorTablePtr_->find(inter.type());

    if (cstrIter == inputsConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown integrator type "
            << inter.type() << nl << nl
            << "Valid integrators are : " << endl
            << inputsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<MultivariateIntegrator<Type>>(cstrIter()(eqn, inter));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::MultivariateIntegrator<Type>::MultivariateIntegrator
(
    const equationType& eqn,
    const dictionary& dict
)
:
    multivariateIntegrator(eqn.nVar(), dict),
    eqnPtr_(&eqn)
{}


template<class Type>
Foam::MultivariateIntegrator<Type>::MultivariateIntegrator
(
    const equationType& eqn,
    const multivariateIntegrator& inter
)
:
    multivariateIntegrator(inter),
    eqnPtr_(&eqn)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::MultivariateIntegrator<Type>::addMidsToInt
(
    const label diri,
    label& fi,
    scalarList& x0,
    scalarList& x1,
    PtrList<Type>& fxs,
    const label li
) const
{
    if (diri >= this->eqnPtr_->nVar())
    {
        scalar dx = x1[0] - x0[0];
        for (label i = 1; i < x0.size(); i++)
        {
            dx *= x1[i] - x0[i];
        }
        fxs.set(fi++, new Type(dx*this->integrateFunc(x0, x1, li)));
        return;
    }

    const scalar xm(0.5*(x1[diri] + x0[diri]));
    const scalar x0Orig(x0[diri]);
    const scalar x1Orig(x1[diri]);

    x1[diri] = xm;
    addMidsToInt(diri+1, fi, x0, x1, fxs, li);

    x0[diri] = xm;
    x1[diri] = x1Orig;
    addMidsToInt(diri+1, fi, x0, x1, fxs, li);
    x0[diri] = x0Orig;

}


template<class Type>
void Foam::MultivariateIntegrator<Type>::integrate_
(
    const PtrList<Type>& Qs,
    const label diri,
    label& fi,
    scalarList& x0,
    scalarList& x1,
    const scalarList& tol,
    const label li,
    Type& fx
) const
{
    if (diri >= this->eqnPtr_->nVar())
    {
        if (!fi)
        {
            fx = integrate_(Qs[fi++], 0, x0, x1, tol, li);
        }
        else
        {
            fx = fx + integrate_(Qs[fi++], 0, x0, x1, tol, li);
        }
        return;
    }

    const scalar xm(0.5*(x1[diri] + x0[diri]));
    const scalar x0Orig(x0[diri]);
    const scalar x1Orig(x1[diri]);

    this->intervals_[diri]++;

    x1[diri] = xm;
    integrate_(Qs, diri+1, fi, x0, x1, tol, li, fx);

    x0[diri] = xm;
    x1[diri] = x1Orig;
    integrate_(Qs, diri+1, fi, x0, x1, tol, li, fx);
    x0[diri] = x0Orig;

}

template<class Type>
Type Foam::MultivariateIntegrator<Type>::integrate_
(
    const Type& Q,
    const label diri,
    const scalarList& X0,
    const scalarList& X1,
    const scalarList& tol,
    const label li
) const
{
    const scalar dx(X1[diri] - X0[diri]);
    if (mag(dx) <= this->minDx_[diri])
    {
        return Q;
    }
    label fi = 0;
    PtrList<Type> fxs(pow(2, X0.size()));
    scalarList x0(X0);
    scalarList x1(X1);
    addMidsToInt
    (
        0,
        fi,
        x0,
        x1,
        fxs,
        li
    );

    Type fx(fxs[0]);
    for (label i = 1; i < fxs.size(); i++)
    {
        fx = fx + fxs[i];
    }
    if (!this->converged(fx, Q, dx, tol[diri], diri))
    {
        fi = 0;
        scalarField eps(tol);
        x0 = X0;
        x1 = X1;

        integrate_
        (
            fxs,
            0,
            fi,
            x0,
            x1,
            eps,
            li,
            fx
        );
    }
    return fx;
}

// ************************************************************************* //
