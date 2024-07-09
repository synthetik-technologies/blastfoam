/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024
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

#include "GaussianQuadrature.H"
#include "polynomialRoots.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::scalar Foam::GaussianQuadrature::binCoeff
(
    const scalar n,
    const scalar k
)
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

void Foam::GaussianQuadrature::calcLegendre
(
    const label n,
    List<scalar>& x,
    List<scalar>& w
)
{
    x.setSize(n);
    w.setSize(n);
    switch (n)
    {
        case 1:
        {
            x[0] = 0.5;
            w[0] = 1.0;
            return;
        }
        case 2:
        {
            x[0] = (1.0 - 1.0/sqrt(3.0))*0.5;
            x[1] = 1.0 - x[0];

            w[0] = 0.5;
            w[1] = 0.5;
            return;
        }
        case 3:
        {
            x[0] = (1.0 - sqrt(3.0/5.0))*0.5;
            x[1] = 0.5;
            x[2] = 1.0 - x[0];

            w[0] = 5.0/18.0;
            w[1] = 8.0/18.0;
            w[2] = w[0];
            return;
        }
        case 4:
        {
            const scalar sqrt65_27(2.0*sqrt(6.0/5.0)/7.0);
            x[0] = (1.0 - sqrt(3.0/7.0 + sqrt65_27))*0.5;
            x[1] = (1.0 - sqrt(3.0/7.0 - sqrt65_27))*0.5;
            x[2] = 1.0 - x[1];
            x[3] = 1.0 - x[0];

            w[0] = (18.0 - sqrt(30.0))/36.0;
            w[1] = (18.0 + sqrt(30.0))/36.0;
            w[2] = w[1];
            w[3] = w[0];
            return;
        }
        case 5:
        {
            const scalar sqrt107_2(2.0*sqrt(10.0/7.0));
            x[0] = (1.0 - sqrt(5.0 + sqrt107_2)/3.0)*0.5;
            x[1] = (1.0 - sqrt(5.0 - sqrt107_2)/3.0)*0.5;
            x[2] = 0.5;
            x[3] = 1.0 - x[1];
            x[4] = 1.0 - x[0];

            w[0] = (322.0 - 13.0*sqrt(70.0))/900.0;
            w[1] = (322.0 + 13.0*sqrt(70.0))/900.0;
            w[2] = 128.0/225.0;
            w[3] = w[1];
            w[4] = w[0];
            return;
        }
        default:
        {
            List<scalar> P(n + 1, 0.0);
            forAll(P, k)
            {
                P[n-k] =
                    binCoeff(n, k)
                   *binCoeff((n + k - 1)/2.0, n)
                   *pow(2.0, n);
            }
            List<scalar> Pp(P.size() - 1);
            forAll(Pp, k)
            {
                Pp[k] = P[k]*scalar(n - k);
            }

            polynomialRoots PRoots(P);
            x = PRoots.rootsRe();
            sort(x);
            forAll(x, i)
            {
                x[i] = (1.0 + x[i])*0.5;
                w[i] =
                    1.0
                   /(
                        (1.0 - sqr(x[i]))
                       *sqr(polynomialRoots::eval(Pp, x[i]))
                    );
            }
        }
    }
}

void Foam::GaussianQuadrature::calcLobatto
(
    const label n,
    List<scalar>& x,
    List<scalar>& w
)
{
    x.setSize(n);
    w.setSize(n);

    switch (n)
    {
        case 2:
        {
            x[0] = 0.0;
            x[1] = 1.0;

            w[0] = 0.5;
            w[1] = 0.5;
            return;
        }
        case 3:
        {
            x[0] = 0.0;
            x[1] = 0.5;
            x[2] = 1.0;

            w[0] = 1.0/6.0;
            w[1] = 4.0/6.0;
            w[2] = 1.0/6.0;
            return;
        }
        case 4:
        {
            x[0] = 0.0;
            x[1] = (1.0 - 0.1*sqrt(5.0))*0.5;
            x[2] = 1.0 - x[1];
            x[3] = 1.0 - x[0];

            w[0] = 1.0/12.0;
            w[1] = 5.0/12.0;
            w[2] = w[1];
            w[3] = w[0];
            return;
        }
        case 5:
        {
            x[0] = 0.0;
            x[1] = (1.0 - sqrt(21.0)/14.0)*0.5;
            x[2] = 0.5;
            x[3] = 1.0 - x[1];
            x[4] = 1.0 - x[0];

            w[0] = 0.05;
            w[1] = 49.0/180.0;
            w[2] = 32.0/90.0;
            w[3] = w[1];
            w[4] = w[0];
            return;
        }
        case 6:
        {
            x[0] = 0.0;
            x[1] = (1.0 - sqrt((7.0 + 2.0*sqrt(7.0))/21.0))*0.5;
            x[2] = (1.0 - sqrt((7.0 - 2.0*sqrt(7.0))/21.0))*0.5;
            x[3] = 1.0 - x[2];
            x[4] = 1.0 - x[1];
            x[5] = 1.0 - x[0];

            w[0] = 1/30.0;
            w[1] = (14.0 - sqrt(7.0))/60.0;
            w[2] = (14.0 + sqrt(7.0))/60.0;
            w[3] = w[2];
            w[4] = w[1];
            w[5] = w[0];
            return;
        }
        default:
        {
            const label n1 = n-1;
            List<scalar> P(n1 + 1, 0.0);
            forAll(P, k)
            {
                P[n1-k] =
                    binCoeff(n1, k)
                   *binCoeff((n1 + k - 1)/2.0, n1)
                   *pow(2.0, n1);
            }

            List<scalar> Pp(P.size() - 1);
            forAll(Pp, k)
            {
                Pp[k] = P[k]*scalar(n1 - k);
            }

            polynomialRoots PpRoots(Pp);
            List<scalar> xPp = PpRoots.rootsRe();
            sort(xPp);

            forAll(xPp, i)
            {
                x[i+1] = (1.0 + xPp[i])*0.5;
                w[i+1] =
                    1.0
                   /(
                        n*(n - 1.0)
                       *sqr(polynomialRoots::eval(P, xPp[i]))
                    );
            }
            x[0] = 0;
            w[0] = 1.0/(n*(n - 1.0));
            x.last() = 1.0;
            w.last() = w[0];
        }
    }
}

// ************************************************************************* //
