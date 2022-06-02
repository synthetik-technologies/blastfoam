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

#include "linearLeastSquares.H"
#include "CoefficientsFwd.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class VarType>
void Foam::linearLeastSquares::setData
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
void Foam::linearLeastSquares::setData
(
    const scalar& x,
    UList<scalar>& xi
)
{
    xi.shallowCopy
    (
        UList<scalar>
        (
            const_cast<scalar*>(&x),
            1
        )
    );
}

template<>
void Foam::linearLeastSquares::setData
(
    const UList<scalar>& x,
    UList<scalar>& xi
)
{
    xi.shallowCopy(x);
}


template<class VarType>
void Foam::linearLeastSquares::findCoeffs
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

    RectangularMatrix<scalar> M(x.size(), eqns.nVar()+1);
    UList<scalar> xi;
    forAll(x, i)
    {
        M(i, 0) = 1.0; // Offset
        setData(x[i], xi);
        forAll(xi, j) // Linear coefficients
        {
            M(i, j+1) = xi[j];
        }
    }
    RectangularMatrix<scalar> MT(M.T());
    coeffs.coeffsRef() = SVDinv(MT*M)*MT*y;
}


template<class VarType>
Foam::autoPtr<Foam::scalarUnivariateEquation> Foam::linearLeastSquares::createEquation
(
    const UList<VarType>& x,
    const scalarField& y,
    const label li
) const
{
    UList<scalar> tmp;
    setData(x.first(), tmp);
    autoPtr<scalarUnivariateEquation> eqns;
    if (tmp.size() == 1)
    {
         eqns.set(new ScalarLinearEquation(tmp.size()));
    }
    else
    {
         eqns.set(new ScalarLinearUnivariateEquation(tmp.size()));
    }
    findCoeffs(eqns(), x, y, li);
    return eqns;
}


template<class VarType>
void Foam::linearLeastSquares::findCoeffs
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

    RectangularMatrix<scalar> M(x.size(), eqns.nVar()+1);
    UList<scalar> xi;
    forAll(x, i)
    {
        M(i, 0) = 1.0; // Offset
        setData(x[i], xi);
        forAll(xi, j) // Linear coefficients
        {
            M(i, j+1) = x[i][j];
        }
    }

    // Apply weights
    RectangularMatrix<scalar> MT(M.T());
    forAll(w, i)
    {
        for (label j=0; j<MT.n(); j++)
        {
            MT(i, j) *= w[i];
        }
    }
    coeffs.coeffsRef() = SVDinv(MT*M)*MT*y;
}

template<class VarType>
Foam::autoPtr<Foam::scalarUnivariateEquation>
Foam::linearLeastSquares::createEquation
(
    const UList<VarType>& x,
    const scalarField& y,
    const scalarField& w,
    const label li
) const
{
    UList<scalar> tmp;
    setData<VarType>(x.first(), tmp);
    autoPtr<scalarUnivariateEquation> eqns;
    if (tmp.size() == 1)
    {
         eqns.set(new ScalarLinearEquation(tmp.size()));
    }
    else
    {
         eqns.set(new ScalarLinearUnivariateEquation(tmp.size()));
    }
    findCoeffs(eqns(), x, y, w, li);
    return eqns;
}

// ************************************************************************* //
