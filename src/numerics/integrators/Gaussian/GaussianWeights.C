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

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

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
        case 6:
        {
            xs =
            {
                -0.932469514203152,
                -0.661209386466265,
                -0.238619186083197,
                0.238619186083197,
                0.661209386466265,
                0.932469514203152
            };
            ws =
            {
                0.171324492379171,
                0.360761573048138,
                0.467913934572690,
                0.467913934572690,
                0.360761573048138,
                0.171324492379171
            };
            break;
        }
        case 7:
        {
            xs =
            {
                -0.949107912342756,
                -0.741531185599397,
                -0.405845151377397,
                0,
                0.405845151377397,
                0.741531185599398,
                0.949107912342754
            };
            ws =
            {
                0.129484966168871,
                0.279705391489287,
                0.381830050505115,
                0.417959183673467,
                0.381830050505115,
                0.279705391489286,
                0.129484966168876
            };
            break;
        }
        case 8:
        {
            xs =
            {
                -0.960289856497539,
                -0.796666477413622,
                -0.525532409916330,
                -0.183434642495650,
                0.183434642495650,
                0.525532409916331,
                0.796666477413621,
                0.960289856497541
            };
            ws =
            {
                0.101228536290378,
                0.222381034453358,
                0.313706645877896,
                0.362683783378362,
                0.362683783378362,
                0.313706645877896,
                0.222381034453358,
                0.101228536290372
            };
            break;
        }
        case 9:
        {
            xs =
            {
                -0.968160239507619,
                -0.836031107326644,
                -0.613371432700588,
                -0.324253423403809,
                0,
                0.324253423403809,
                0.613371432700589,
                0.836031107326646,
                0.968160239507618
            };
            ws =
            {
                8.127438836158031e-02,
                1.806481606948867e-01,
                2.606106964029209e-01,
                3.123470770400045e-01,
                3.302393550012606e-01,
                3.123470770400047e-01,
                2.606106964029202e-01,
                1.806481606948852e-01,
                8.127438836158338e-02
            };
            break;
        }
        case 10:
        {
            xs =
            {
                -0.973906528517197,
                -0.865063366688947,
                -0.679409568299039,
                -0.433395394129245,
                -0.148874338981631,
                0.148874338981631,
                0.433395394129245,
                0.679409568299040,
                0.865063366688948,
                0.973906528517198
            };
            ws =
            {
                6.667134430865264e-02,
                1.494513491505223e-01,
                2.190863625160392e-01,
                2.692667193099844e-01,
                2.955242247147547e-01,
                2.955242247147547e-01,
                2.692667193099837e-01,
                2.190863625160358e-01,
                1.494513491505162e-01,
                6.667134430865065e-02
            };
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unsupported number of Nodes for Gaussian integration" << nl
                << "1 to 10 nodes can be used."
                << abort(FatalError);
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
