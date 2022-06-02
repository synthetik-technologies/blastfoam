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

#include "nonLinearLeastSquares.H"
#include "CoefficientsFwd.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonLinearLeastSquares, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonLinearLeastSquares::nonLinearLeastSquares
(
    const scalar tol,
    const label maxSteps
)
:
    tolerance_(tol),
    maxSteps_(maxSteps)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nonLinearLeastSquares::findCoeffs
(
    scalarUnivariateEquation& eqns,
    const List<scalarList>& x,
    const scalarList& y,
    const label li
) const
{
    scalarCoefficients& coeffs
    (
        dynamicCast<scalarCoefficients>(eqns)
    );
    RectangularMatrix<scalar> J(x.size(), coeffs.n());
    RectangularMatrix<scalar> JT(J.T());
    scalarField r(x.size());
    scalar eOld = great;

    for (label stepi = 0; stepi < maxSteps_; stepi++)
    {
        scalar e = 0.0;
        forAll(r, i)
        {
            r[i] = y[i] - eqns.fX(x[i], li);
            e += sqr(r[i]);
        }
        e = sqrt(e)/scalar(x.size());
        coeffs.coeffJ(x, li, J);
        JT = J.T();
        coeffs.coeffsRef() -= SVDinv(JT*J)*JT*r;

        if (mag(e - eOld) < tolerance_)
        {
            break;
        }
        eOld = e;
    }
}


void Foam::nonLinearLeastSquares::findCoeffs
(
    scalarUnivariateEquation& eqns,
    const List<scalarList>& x,
    const scalarList& y,
    const scalarList& w,
    const label li
) const
{
    scalarCoefficients& coeffs
    (
        dynamicCast<scalarCoefficients>(eqns)
    );
    RectangularMatrix<scalar> J(x.size(), coeffs.n());
    RectangularMatrix<scalar> JT(J.T());
    scalarField r(x.size());
    scalar eOld = great;

    for (label stepi = 0; stepi < maxSteps_; stepi++)
    {
        scalar e = 0.0;
        forAll(r, i)
        {
            r[i] = y[i] - eqns.fX(x[i], li);
            e += sqr(r[i]);
        }
        e = sqrt(e)/scalar(x.size());
        coeffs.coeffJ(x, li, J);
        JT = J.T();
        forAll(w, i)
        {
            for (label j=0; j<JT.n(); j++)
            {
                JT(i, j) *= w[i];
            }
        }
        coeffs.coeffsRef() -= SVDinv(JT*J)*JT*r;

        if (mag(e - eOld) < tolerance_)
        {
            break;
        }
        eOld = e;
    }
}

// ************************************************************************* //
