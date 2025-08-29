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

#include "AntoineExtendedSaturationPressureModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationPressureModels
{
    defineTypeNameAndDebug(AntoineExtended, 0);
    addToRunTimeSelectionTable
    (
        saturationPressureModel,
        AntoineExtended,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationPressureModels::AntoineExtended::AntoineExtended
(
    const dictionary& dict
)
:
    saturationPressureModel(dict),
    A_("A", dimless, dict),
    B_("B", dimTemperature, dict),
    C_("C", dimTemperature, dict),
    D_("D", inv(dimTemperature), dict),
    E_("E", inv(pow3(dimTemperature)), dict),
    F_("F", inv(sqr(dimTemperature)), dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationPressureModels::AntoineExtended::~AntoineExtended()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::saturationPressureModels::AntoineExtended::pSat
(
    const scalar T
) const
{
    const scalar TC = max(T - Toffset_.value(), 0.0) + C_.value();
    return
        pow
        (
            10.0,
            A_.value()
          - B_.value()/TC
          + D_.value()*TC
          + E_.value()*sqr(TC)
          + F_.value()*pow3(TC)
        );
}


Foam::scalar Foam::saturationPressureModels::AntoineExtended::derivative
(
    const scalar T
) const
{
    static const scalar ln10 = log(10.0);
    const scalar TC = max(T - Toffset_.value(), 0.0) + C_.value();
    return
        (
            B_.value()/sqr(TC)
          + D_.value()
          + 2.0*E_.value()*TC
          + 3.0*F_.value()*sqr(TC)
        )
       *ln10
       *pow
        (
            10.0,
            A_.value()
          - B_.value()/TC
          + D_.value()*TC
          + E_.value()*sqr(TC)
          + F_.value()*pow3(TC)
        );
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::saturationPressureModels::AntoineExtended::pSat
(
    const volScalarField::Internal& T
) const
{
    const volScalarField::Internal TC(max(T - Toffset_, zeroT) + C_);
    return
        pow
        (
            10.0,
            A_ - B_/TC + TC*(D_ + TC*(E_ + TC*F_))
        );
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::saturationPressureModels::AntoineExtended::derivative
(
    const volScalarField::Internal& T
) const
{
    static const scalar ln10 = log(10.0);
    const volScalarField::Internal TC(max(T - Toffset_, zeroT) + C_);

    return
        (B_/sqr(TC) + D_ + TC*(2.0*E_ + 3.0*TC*F_))
       *ln10
       *pow(10.0, A_ - B_/TC + TC*(D_ + TC*(E_ + TC*F_)));
}


// ************************************************************************* //
