/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "polynomialRoots.H"
#include "eigenSolver.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polynomialRoots::polynomialRoots
(
    const List<scalar>& coeffs
)
:
    coeffs_(coeffs),
    rootsRe_(coeffs_.size() - 1),
    rootsIm_(coeffs_.size() - 1)
{
    scalarSquareMatrix A(coeffs_.size()-1, 0.0);
    scalar c0 = stabilise(coeffs_[0], small);
    for (label i = 1; i < coeffs_.size(); i++)
    {
        A(0, i-1) = -coeffs_[i]/c0;
    }
    for (label i = 1; i < coeffs_.size()-1; i++)
    {
        A(i, i-1) =  1.0;
    }
    eigenSolver eig(A, false);
    rootsRe_ = eig.eigenvaluesRe();
    rootsIm_ = eig.eigenvaluesIm();

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::polynomialRoots::eval(const scalar x) const
{
    return eval(coeffs_, x);
}


Foam::scalar Foam::polynomialRoots::eval
(
    const List<scalar>& coeffs,
    const scalar x
)
{
    label i = coeffs.size()-1;
    scalar res = coeffs[i--];
    scalar xPow = x;

    res += xPow*coeffs[i--];
    while (i >= 0)
    {
        xPow *= x;
        res += xPow*coeffs[i--];
    }
    return res;
}


Foam::word Foam::polynomialRoots::polyName(const List<scalar>& coeffs)
{
    label i = coeffs.size()-1;
    OStringStream os;
    scalar c = coeffs[i--];
    if (mag(c) > 0)
    {
        os << c;
    }
    else if (i < 0)
    {
        return "0";
    }
    if (i < 0)
    {
        return os.str();
    }

    c = coeffs[i--];
    if (c > 0)
    {
        if (c == 1)
        {
            os << " + x";
        }
        else
        {
            os << " + " << c << "*x";
        }
    }
    else if (c < 0)
    {
        if (c == -1)
        {
            os << " + x";
        }
        else
        {
           os << " - " << mag(c)<< "*x";
        }
    }
    if (i < 0)
    {
        return os.str();
    }

    label p = 2;
    while (i >= 0)
    {
        c = coeffs[i--];
        if (c > 0)
        {
            if (c == 1)
            {
                os << " + pow(x," << p++ << ")";
            }
            else
            {
                os << " + " << c << "*pow(x," << p++ << ")";
            }

        }
        else if (c < 0)
        {
            if (c == -1)
            {
                os << " + pow(x," << p++ << ")";
            }
            else
            {
              os << " - " << mag(c)<< "*pow(x," << p++ << ")";
            }
        }
    }
    return os.str();
}


Foam::word Foam::polynomialRoots::polyName() const
{
    return polyName(coeffs_);
}


Foam::List<Foam::scalar> Foam::polynomialRoots::derivative
(
    const List<scalar>& P
)
{
    if (P.size() == 1)
    {
        return List<scalar>(1, 0.0);
    }
    
    List<scalar> Pprime(P.size() - 1);
    label n = P.size() - 1;
    forAll(Pprime, k)
    {
        Pprime[k] = P[k]*scalar(n-k);
    }
    return Pprime;
}


Foam::List<Foam::scalar> Foam::polynomialRoots::derivative() const 
{
    return derivative(coeffs_);
}


// ************************************************************************* //
