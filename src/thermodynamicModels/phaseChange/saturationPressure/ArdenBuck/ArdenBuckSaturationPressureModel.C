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

#include "ArdenBuckSaturationPressureModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationPressureModels
{
    defineTypeNameAndDebug(ArdenBuck, 0);
    addToRunTimeSelectionTable(saturationPressureModel, ArdenBuck, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationPressureModels::ArdenBuck::ArdenBuck(const dictionary& dict)
:
    saturationPressureModel(dict),
    A_("A", dimPressure, dict.lookupOrDefault("A", 611.21)),
    B_("B", dimless, dict.lookupOrDefault("B", 18.678)),
    C_("C", dimTemperature, dict.lookupOrDefault("C", 234.5)),
    D_("D", dimTemperature, dict.lookupOrDefault("D", 257.14))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationPressureModels::ArdenBuck::~ArdenBuck()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::saturationPressureModels::ArdenBuck::pSat
(
    const scalar T
) const
{
    const scalar Toff = max(T - Toffset_.value(), 0.0);
    return
        A_.value()
       *exp((B_.value() - Toff/C_.value())*Toff/(D_.value() + Toff));
}


Foam::scalar Foam::saturationPressureModels::ArdenBuck::derivative
(
    const scalar T
) const
{
    const scalar Toff = max(T - Toffset_.value(), 0.0);
    const scalar x = (B_.value() - Toff)/(C_.value()*(D_.value() + Toff));
    const scalar dxdT =
      - (B_.value()*C_.value() + D_.value())
       /(C_.value()*sqr(D_.value() + Toff));
    return A_.value()*exp(Toff*x)*(x - Toff*dxdT);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::saturationPressureModels::ArdenBuck::pSat
(
    const volScalarField::Internal& T
) const
{
    const volScalarField::Internal Toff(max(T - Toffset_, zeroT));
    return A_*exp((B_ - Toff/C_)*Toff/(D_ + Toff));
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::saturationPressureModels::ArdenBuck::derivative
(
    const volScalarField::Internal& T
) const
{
    const volScalarField::Internal Toff(max(T - Toffset_, zeroT));
    const volScalarField::Internal x((B_ - Toff)/(C_*(D_ + Toff)));
    return
        A_*exp(Toff*x)
       *(x - Toff*((B_*C_ + D_)/(C_*sqr(D_ + Toff))));
}

// ************************************************************************* //
