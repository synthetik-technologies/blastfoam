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

#include "AntoineSaturationPressureModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationPressureModels
{
    defineTypeNameAndDebug(Antoine, 0);
    addToRunTimeSelectionTable(saturationPressureModel, Antoine, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationPressureModels::Antoine::Antoine(const dictionary& dict)
:
    saturationPressureModel(dict),
    A_("A", dimless, dict.lookup("A")),
    B_("B", dimTemperature, dict.lookup("B")),
    C_("C", dimTemperature, dict.lookup("C"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationPressureModels::Antoine::~Antoine()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::saturationPressureModels::Antoine::pSat
(
    const scalar T
) const
{
    const scalar Toff = max(T - Toffset_.value(), 0.0);
    return
        pow(10.0, A_.value() - B_.value()/(Toff + C_.value()));
}


Foam::scalar Foam::saturationPressureModels::Antoine::derivative
(
    const scalar T
) const
{
    static const scalar ln10 = log(10.0);
    const scalar Toff = max(T - Toffset_.value(), 0.0);
    return
        ln10*B_.value()/sqr(Toff + C_.value())
       *pow(10.0, A_.value() - B_.value()/(Toff + C_.value()));
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::saturationPressureModels::Antoine::pSat
(
    const volScalarField::Internal& T
) const
{
    return pow(10.0, A_ - B_/(max(T - Toffset_, zeroT) + C_));
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::saturationPressureModels::Antoine::derivative
(
    const volScalarField::Internal& T
) const
{
    static const scalar ln10 = log(10.0);
    volScalarField::Internal TC(max(T - Toffset_, zeroT) + C_);

    return ln10*B_/sqr(TC)*pow(10.0, A_ - B_/TC);
}

// ************************************************************************* //
