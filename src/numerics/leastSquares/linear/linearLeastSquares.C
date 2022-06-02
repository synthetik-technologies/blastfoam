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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearLeastSquares, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::linearLeastSquares::linearLeastSquares()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::linearLeastSquares::findCoeffs
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

    RectangularMatrix<scalar> M(x.size(), eqns.nVar()+1);
    forAll(x, i)
    {
        M(i, 0) = 1.0; // Offset
        forAll(x[i], j) // Linear coefficients
        {
            M(i, j+1) = x[i][j];
        }
    }
    RectangularMatrix<scalar> MT(M.T());
    coeffs.coeffsRef() = SVDinv(MT*M)*MT*scalarField(y);
}


void Foam::linearLeastSquares::findCoeffs
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

    RectangularMatrix<scalar> M(x.size(), eqns.nVar()+1);
    forAll(x, i)
    {
        M(i, 0) = 1.0; // Offset
        forAll(x[i], j) // Linear coefficients
        {
            M(i, j+1) = x[i][j];
        }
    }
    RectangularMatrix<scalar> MT(M.T());
    forAll(w, i)
    {
        for (label j=0; j<MT.n(); j++)
        {
            MT(i, j) *= w[i];
        }
    }
    coeffs.coeffsRef() = SVDinv(MT*M)*MT*scalarField(y);
}



// ************************************************************************* //
