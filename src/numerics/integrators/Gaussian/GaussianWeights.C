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

#include "GaussianWeights.H"
#include "ListOps.H"
#include "polynomialRoots.H"

// * * * * * * * * * * * * Static Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::GaussianWeights::binCoeff(const scalar n, const scalar k)
{
    scalar nmk(n - k + 1);
    if (nmk < small)
    {

        label nmki = round(nmk);
        if (mag(nmk - nmki) < small)
        {
            return 0.0;
        }
    }
    return tgamma(n + 1.0)/(tgamma(k + 1.0)*tgamma(nmk));
}


void Foam::GaussianWeights::calcQuadrature
(
    const label N,
    List<scalar>& ws,
    List<scalar>& xs
)
{
    ws.resize(N);
    xs.resize(N);

    List<scalar> P(N+1, 0.0);
    forAll(P, k)
    {
        P[N-k] = binCoeff(N, k)*binCoeff((N+k-1)/2.0, N)*pow(2.0, N);
    }
    List<scalar> Pprime(P.size() - 1);
    forAll(Pprime, k)
    {
        Pprime[k] = P[k]*scalar(N-k);
    }

    polynomialRoots PRoots(P);
    xs = PRoots.rootsRe();
    sort(xs);
    forAll(xs, i)
    {
        ws[i] =
            2.0
           /(
                (1.0 - sqr(xs[i]))
               *sqr(polynomialRoots::eval(Pprime, xs[i]))
            );
    }
}


void Foam::GaussianWeights::setQuadrature
(
    const label N,
    List<scalar>& ws,
    List<scalar>& xs
)
{
    switch (N)
    {
        case 1:
        {
            xs = {0.0};
            ws = {2.0};
            break;
        }
        case 2:
        {
            const scalar x = 1.0/sqrt(3.0);
            xs = {-x, x};
            ws = {1.0, 1.0};
            break;
        }
        case 3:
        {
            const scalar x = sqrt(3.0/5.0);
            xs = {-x, 0.0,  x};
            ws = {5.0/9.0, 8.0/9.0, 5.0/9.0};
            break;
        }
        case 4:
        {
            const scalar x1 = sqrt(3.0/7.0 - 2.0/7.0*sqrt(6.0/5.0));
            const scalar x2 = sqrt(3.0/7.0 + 2.0/7.0*sqrt(6.0/5.0));

            const scalar w1 = (18.0 + sqrt(30.0))/36.0;
            const scalar w2 = (18.0 - sqrt(30.0))/36.0;

            xs = {-x2, -x1, x1, x2};
            ws = {w2, w1, w1, w2};
            break;
        }
        case 5:
        {
            const scalar x1 = sqrt(5.0 - 2.0*sqrt(10.0/7.0))/3.0;
            const scalar x2 = sqrt(5.0 + 2.0*sqrt(10.0/7.0))/3.0;

            const scalar w0 = 128.0/225.0;
            const scalar w1 = (322.0 + 13.0*sqrt(70.0))/900.0;
            const scalar w2 = (322.0 - 13.0*sqrt(70.0))/900.0;

            xs = {-x2, -x1, 0.0, x1, x2};
            ws = {w2, w1, w0, w1, w2};
            break;
        }
        default:
        {
            calcQuadrature(N, ws, xs);
            break;
        }
    }

    // Normalize weights and move xs to be in (0, 1)
    scalar sumW = 0.0;
    forAll(ws, i)
    {
        sumW += ws[i];
    }
    forAll(xs, i)
    {
        ws[i] /= sumW;
        xs[i] = 0.5*(xs[i] + 1.0);
    }
}

// ************************************************************************* //
