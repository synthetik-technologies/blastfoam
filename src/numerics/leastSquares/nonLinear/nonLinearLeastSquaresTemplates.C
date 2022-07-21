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

#include "nonLinearLeastSquares.H"
#include "coefficientsFwd.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

namespace Foam
{
template<class VarType>
void nonLinearLeastSquares::setData
(
    const VarType& x,
    UList<scalar>& xi
)
{
    xi.shallowCopy
    (
        UList<scalar>
        (
            const_cast<scalar*>(x.v_),
            pTraits<VarType>::nComponents
        )
    );
}


template<>
void nonLinearLeastSquares::setData
(
    const scalar& x,
    UList<scalar>& xi
)
{
    xi.shallowCopy
    (
        UList<scalar>(const_cast<scalar*>(&x), 1)
    );
}

template<>
void nonLinearLeastSquares::setData
(
    const UList<scalar>& x,
    UList<scalar>& xi
)
{
    xi.shallowCopy(x);
}

} // End namespace Foam


template<class VarType>
void Foam::nonLinearLeastSquares::findCoeffs
(
    scalarUnivariateEquation& eqns,
    const UList<VarType>& x,
    const scalarField& y,
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

    UList<scalar> xi;
    for (label stepi = 0; stepi < maxSteps_; stepi++)
    {
        scalar e = 0.0;
        forAll(r, i)
        {
            setData(x[i], xi);

            r[i] = y[i] - eqns.fX(xi, li);
            e += sqr(r[i]);

            coeffs.coeffJ(xi, i, J, li);
        }
        e = sqrt(e)/scalar(x.size());

        JT = J.T();
        coeffs.coeffsRef() -= SVDinv(JT*J)*JT*r;

        if (mag(e - eOld) < tolerance_)
        {
            break;
        }
        eOld = e;
    }
}


template<class VarType>
void Foam::nonLinearLeastSquares::findCoeffs
(
    scalarUnivariateEquation& eqns,
    const UList<VarType>& x,
    const scalarField& y,
    const scalarField& w,
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

    UList<scalar> xi;
    for (label stepi = 0; stepi < maxSteps_; stepi++)
    {
        scalar e = 0.0;
        forAll(r, i)
        {
            setData(x[i], xi);

            r[i] = y[i] - eqns.fX(xi, li);
            e += sqr(r[i]);

            coeffs.coeffJ(xi, i, J, li);
        }
        e = sqrt(e)/scalar(x.size());

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
