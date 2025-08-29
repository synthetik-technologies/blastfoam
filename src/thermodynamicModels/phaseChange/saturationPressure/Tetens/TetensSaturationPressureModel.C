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

#include "TetensSaturationPressureModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationPressureModels
{
    defineTypeNameAndDebug(Tetens, 0);
    addToRunTimeSelectionTable(saturationPressureModel, Tetens, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationPressureModels::Tetens::Tetens(const dictionary& dict)
:
    saturationPressureModel(dict),
    A_("A", dimPressure, dict.lookupOrDefault("A", 610.78)),
    B_("B", dimless, dict.lookupOrDefault("B", 17.27)),
    C_("C", dimTemperature, dict.lookupOrDefault("C", 243.04))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationPressureModels::Tetens::~Tetens()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::saturationPressureModels::Tetens::pSat
(
    const scalar T
) const
{
    const scalar Toff = max(T - Toffset_.value(), 0.0);
    return A_.value()*exp(B_.value()*Toff/(Toff + C_.value()));
}


Foam::scalar Foam::saturationPressureModels::Tetens::derivative
(
    const scalar T
) const
{
    const scalar Toff = max(T - Toffset_.value(), 0.0);
    return
        exp(B_.value()*Toff/(Toff + C_.value()))
       *A_.value()*B_.value()*C_.value()
       /sqr(Toff + C_.value());
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::saturationPressureModels::Tetens::pSat
(
    const volScalarField::Internal& T
) const
{
    const volScalarField::Internal Toff(max(T - Toffset_, zeroT));
    return A_*exp(B_*Toff/(Toff + C_));
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::saturationPressureModels::Tetens::derivative
(
    const volScalarField::Internal& T
) const
{
    const volScalarField::Internal Toff(max(T - Toffset_, zeroT));
    return exp(B_*Toff/(Toff + C_))*A_*B_*C_/sqr(Toff + C_);
}

// ************************************************************************* //
