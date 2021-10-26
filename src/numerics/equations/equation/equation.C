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

#include "equation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::equation::equation
(
    const scalar lowerLimit,
    const scalar upperLimit
)
:
    Equation<scalar, scalar>(1, 1, lowerLimit, upperLimit)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::equation::~equation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::equation::nVar() const
{
    return Equation<scalar, scalar>::nVar();
}


Foam::label Foam::equation::nEqns() const
{
    return Equation<scalar, scalar>::nEqns();
}


Foam::tmp<Foam::scalarField> Foam::equation::lowerLimits() const
{
    return tmp<scalarField>(new scalarField(1, this->lowerLimits_));
}


Foam::tmp<Foam::scalarField> Foam::equation::upperLimits() const
{
    return tmp<scalarField>(new scalarField(1, this->upperLimits_));
}


Foam::tmp<Foam::scalarField> Foam::equation::dX() const
{
    return tmp<scalarField>(new scalarField(1, this->dx_));
}


void Foam::equation::setDX(const scalarField& newDx) const
{
    this->dx_ = newDx[0];
}


void Foam::equation::limit(scalar& x) const
{
    Equation<scalar, scalar>::limit(x);
}


void Foam::equation::f
(
    const scalar& x,
    const label li,
    scalar& fx
) const
{
    fx = this->fx(x, li);
}


void Foam::equation::f
(
    const scalarField& x,
    const label li,
    scalar& fx
) const
{
    fx = this->fx(x[0], li);
}


void Foam::equation::f
(
    const scalarField& x,
    const label li,
    scalarField& fx
) const
{
    fx[0] = this->fx(x[0], li);
}


void Foam::equation::gradient
(
    const scalarField& x,
    const label li,
    scalar& fx,
    scalarField& grad
) const
{
    fx = this->fx(x[0], li);
    grad[0] = this->dfdx(x[0], li);
}


void Foam::equation::jacobian
(
    const scalarField& x,
    const label li,
    scalarField& fx,
    RectangularMatrix<scalar>& J
) const
{
    fx[0] = this->fx(x[0], li);
    J(0, 0) = this->dfdx(x[0], li);
}


bool Foam::equation::containsRoot
(
    const scalar y0,
    const scalar y1
) const
{
    if (y0*y1 > 0)
    {
        #ifdef FULLDEBUG
        FatalErrorInFunction
            << "Solution is not bracked:" << nl
            << "limits: (" << lower() << ","<< upper() << ")" << endl
            << "f(x0)=" << y0 << ", f(x1)=" << y1 << endl
            << abort(FatalError);
        #endif
        return false;
    }
    return true;
}


bool Foam::equation::containsRoot
(
    const scalarField& y0s,
    const scalarField& y1s
) const
{
    if (y0s[0]*y1s[0] > 0)
    {
        #ifdef FULLDEBUG
        FatalErrorInFunction
            << "Solution is not bracked:" << nl
            << "limits: (" << lower() << ","<< upper() << ")" << endl
            << "f(x0)=" << y0s[0] << ", f(x1)=" << y1s[0] << endl
            << abort(FatalError);
        #endif
        return false;
    }
    return true;
}


bool Foam::equation::containsRoot(const label li) const
{
    return containsRoot
    (
        this->fx(lower(), li),
        this->fx(upper(), li)
    );
}


// ************************************************************************* //
