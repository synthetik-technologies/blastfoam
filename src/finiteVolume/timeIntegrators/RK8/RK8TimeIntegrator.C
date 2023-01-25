/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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

#include "RK8TimeIntegrator.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace timeIntegrators
{
    defineTypeNameAndDebug(RK8, 0);
    addToRunTimeSelectionTable(timeIntegrator, RK8, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrators::RK8::RK8
(
    const fvMesh& mesh,
    Istream& is
)
:
    timeIntegrator(mesh)
{
    this->as_ =
    {
        {1.0},
        {1.0, 0.0},
        {1.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
    };

    const scalar sqrt21 = sqrt(21.0);
    const scalar b_12_8 = 49.0/180.0;
    const scalar b_12_7 = -b_12_8 + 49.0/180;

    const scalar b_10_5 = 1.0/9.0;

    scalar b_12_1 = 1.0/20.0;
    scalar b_12_2 = 0.0;
    scalar b_12_3 = 0.0;
    scalar b_12_4 = 0.0;
    scalar b_12_5 = 0.0;
    scalar b_12_6 = 0.0;
    scalar b_12_9 = 16.0/45;
    scalar b_12_10 = 49.0/180.0;
    scalar b_12_11 = 1.0/20.0;

    scalar c_2 = 0.5;
    scalar c_3 = 0.5;
    scalar c_4 = (7.0 + sqrt21)/14.0;
    scalar c_5 = (7.0 + sqrt21)/14.0;
    scalar c_6 = 0.5;
    scalar c_7 = (7.0 - sqrt21)/14.0;
    scalar c_8 = (7.0 - sqrt21)/14.0;
    // scalar c_9 = 0.5;
    // scalar c_10 = (7.0 + sqrt21)/14.0;
    // scalar c_11 = 1.0;

    scalar b_2_1 = 0.5;

    scalar b_3_1 = 0.25;
    scalar b_3_2 = 0.25;

    scalar b_4_1 = 1.0/7.0;
    scalar b_4_2 = (-7.0 - 3.0*sqrt21)/98.0;
    scalar b_4_3 = (21.0 + 5.0*sqrt21)/49.0;

    scalar b_5_1 = (11.0 + sqrt21)/84.0;
    scalar b_5_2 = 0.0;
    scalar b_5_3 = 4.0*sqrt21/63.0 + 2.0/7.0;
    scalar b_5_4 = (21.0 - sqrt21)/252.0;

    scalar b_6_1 = (5.0 + sqrt21)/48.0;
    scalar b_6_2 = 0.0;
    scalar b_6_3 = (9.0 + sqrt21)/36.0;
    scalar b_6_4 = (-231.0 + 14.0*sqrt21)/360.0;
    scalar b_6_5 = (63.0 - 7.0*sqrt21)/80.0;

    scalar b_7_1 = (10.0 - sqrt21)/42.0;
    scalar b_7_2 = 0.0;
    scalar b_7_3 =
      - 24.0/35.0*b_10_5 - 136.0/105.0 - 12.0/245.0*b_10_5*sqrt21 + 656.0/2205.0*sqrt21;
    scalar b_7_4 =
        7.0 - 3.0/10.0*b_10_5*sqrt21 - 71.0/45.0*sqrt21 + 3.0/10.0*b_10_5;
    scalar b_7_5 =
        -3.0/10.0*b_10_5 + 3.0/10.0*b_10_5*sqrt21 - 43.0/6.0 + 169.0/105.0*sqrt21;
    scalar b_7_6 =
      - 277.0/735.0*sqrt21 + 181.0/105.0 + 12.0/245.0*b_10_5*sqrt21 + 24.0/35.0*b_10_5;

    scalar b_8_1 =
        -(180*b_12_8*sqrt21 - 49.0*sqrt21 - 1800.0*b_12_8 + 343.0)/7560.0*b_12_8;
    scalar b_8_2 = 0.0;
    scalar b_8_5 =
      - (
            441.0*b_10_5*sqrt21
          - 3240.0*b_7_5*b_12_8
          - 28.0*sqrt21
          + 882.0*b_7_5
          - 2205.0*b_10_5
          + 147.0
        )/(3240.0*b_12_8);
    scalar b_8_6 =
        (
            72.0*b_10_5*sqrt21
          + 1620.0*b_7_6*b_12_8
          - 29.0*sqrt21
          - 441.0*b_7_6
          - 252.0*b_10_5
          + 119.0
        )/(1620.0*b_12_8);
    scalar b_8_3 =
      - (
            900.0*b_12_8*sqrt21
          + 11340.0*b_7_2*b_12_8
          + 11340.0*b_8_6*b_12_8
          - 98.0*sqrt21
          - 3087.0*b_7_2
          - 4860.0*b_12_8
          + 686.0
        )/(11340.0*b_12_8);

    scalar b_8_7 = 49.0/(1620.0*b_12_8);
    scalar b_8_4 =
        (sqr(c_8)/2.0 - b_8_2*c_2 - b_8_3*c_3 - b_8_5*c_5 - b_8_6*c_6 - b_8_7*c_7)/c_4;

    scalar b_9_1 = 1.0/32.0;
    scalar b_9_2 = 0.0;
    scalar b_9_3 =
        1.0/8.0*b_10_5*sqrt21 - 1.0/8.0*b_10_5 - 1.0/72.0*sqrt21 + 1.0/72.0;
    scalar b_9_4 =
      - 49.0/288.0 - 7.0/32.0*b_10_5*sqrt21 + 7.0/288.0*sqrt21 + 49.0/32.0*b_10_5;
    scalar b_9_5 =
        7.0/32.0*b_10_5*sqrt21 - 35.0/576.0*sqrt21 - 49.0/32.0*b_10_5 + 21.0/64.0;
    scalar b_9_6 =
      - 1.0/8.0*b_10_5*sqrt21 + 1.0/8.0*b_10_5 + 1.0/72.0*sqrt21 + 5.0/36.0;
    scalar b_9_7 =
        91.0/576.0 + 7.0/192.0*sqrt21 - 585.0/1568.0*b_12_8*sqrt21 - 405.0/224.0*b_12_8;
    scalar b_9_8 = 585.0/1568.0*sqrt21*b_12_8 + 405.0/224.0*b_12_8;

    scalar b_10_1 = 1.0/14.0;
    scalar b_10_2 = 0.0;
    scalar b_10_9 = 4.0*sqrt21/35.0 + 132.0/245.0;
    scalar b_10_3 =
        -6.0/49.0*b_10_5*sqrt21 - 2.0/7.0*b_10_5 + 2.0/147.0*sqrt21 + 2.0/63.0;
    scalar b_10_4 = 1.0/9.0 - b_10_5;
    scalar b_10_6 =
        2.0/7.0*b_10_5 - 803.0/2205.0 + 6.0/49.0*b_10_5*sqrt21 - 59.0/735.0*sqrt21;
    scalar b_10_7 =
        1.0/9.0 + 1.0/42.0*sqrt21 + 2295.0/686.0*b_12_8 + 495.0/686.0*b_12_8*sqrt21;
    scalar b_10_8 = (-2295.0/686.0*b_12_8 - 495.0/686.0)*b_12_8*sqrt21;

    scalar b_11_1 = 0.0;
    scalar b_11_2 = 0.0;
    scalar b_11_9 = (28.0 - 28*sqrt21)/45.0;
    scalar b_11_10 = (49.0 - 7.0*sqrt21)/18.0;
    scalar b_11_3 =
        2.0/3.0*b_10_5*sqrt21 - 2.0/3.0*b_10_5 - 2.0/27.0*sqrt21 + 2.0/27.0;
    scalar b_11_4 =
        -7.0/6.0*b_10_5*sqrt21 + 7.0/54.0*sqrt21 + 49.0/6.0*b_10_5 - 49.0/54.0;
    scalar b_11_5 =
        7.0/27.0*sqrt21 - 77.0/54.0 - 49.0/6.0*b_10_5 + 7.0/6.0*b_10_5*sqrt21;
    scalar b_11_6 =
        2.0/3.0*b_10_5 - 64.0/135.0 - 2.0/3.0*b_10_5*sqrt21 + 94.0/135.0*sqrt21;
    scalar b_11_7 =
        7.0/18.0 - 265.0/98.0*b_12_8*sqrt21 - 215.0/14.0*b_12_8;
    scalar b_11_8 = 265.0/98.0*b_12_8*sqrt21 + 215.0/14.0*b_12_8;

    this->bs_ =
    {
        {b_2_1},
        {b_3_1, b_3_2},
        {b_4_1, b_4_2, b_4_3},
        {b_5_1, b_5_2, b_5_3, b_5_4},
        {b_6_1, b_6_2, b_6_3, b_6_4, b_6_5},
        {b_7_1, b_7_2, b_7_3, b_7_4, b_7_5, b_7_6},
        {b_8_1, b_8_2, b_8_3, b_8_4, b_8_5, b_8_6, b_8_7},
        {b_9_1, b_9_2, b_9_3, b_9_4, b_9_5, b_9_6, b_9_7, b_9_8},
        {b_10_1, b_10_2, b_10_3, b_10_4, b_10_5, b_10_6, b_10_7, b_10_8, b_10_9},
        {b_11_1, b_11_2, b_11_3, b_11_4, b_11_5, b_11_6, b_11_7, b_11_8, b_11_9, b_11_10},
        {b_12_1, b_12_2, b_12_3, b_12_4, b_12_5, b_12_6, b_12_7, b_12_8, b_12_9, b_12_10, b_12_11}
    };
    this->initialize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrators::RK8::~RK8()
{}
// ************************************************************************* //
