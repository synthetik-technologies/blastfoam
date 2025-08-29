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

#include "constantSaturationPressureModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationPressureModels
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable(saturationPressureModel, constant, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationPressureModels::constant::constant(const dictionary& dict)
:
    saturationPressureModel(dict),
    pSat_("pSat", dimPressure, dict.lookup("pSat"))
{}


Foam::saturationPressureModels::constant::constant
(
    const dictionary& dict,
    const scalar p
)
:
    saturationPressureModel(dict),
    pSat_("pSat", dimPressure, p)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationPressureModels::constant::~constant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::saturationPressureModels::constant::pSat
(
    const scalar T
) const
{
    return pSat_.value();
}


Foam::scalar Foam::saturationPressureModels::constant::derivative
(
    const scalar T
) const
{
    return 0.0;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::saturationPressureModels::constant::pSat
(
    const volScalarField::Internal& T
) const
{
    return volScalarField::Internal::New
    (
        IOobject::groupName("pSat", T.group()),
        T.mesh(),
        pSat_
    );
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::saturationPressureModels::constant::derivative
(
    const volScalarField::Internal& T
) const
{
    return volScalarField::Internal::New
    (
        IOobject::groupName("dpSatdT", T.group()),
        T.mesh(),
        dimensionedScalar(dimPressure/dimTemperature, 0.0)
    );
}

// ************************************************************************* //
