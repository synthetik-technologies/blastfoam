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

#include "functionSaturationPressureModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationPressureModels
{
    defineTypeNameAndDebug(function, 0);
    addToRunTimeSelectionTable(saturationPressureModel, function, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationPressureModels::function::function(const dictionary& dict)
:
    saturationPressureModel(dict),
    func_
    (
        Function1<scalar>::New
        (
            "pSatFunction",
            dimTemperature,
            dimPressure,
            dict
        )
    )
{
    if (!dict.found("Toffset"))
    {
        Toffset_.value() = 0.0;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationPressureModels::function::~function()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::saturationPressureModels::function::pSat
(
    const scalar T
) const
{
    return func_->value(T - Toffset_.value());
}


Foam::scalar Foam::saturationPressureModels::function::derivative
(
    const scalar T
) const
{
    const scalar T0 = T - Toffset_.value();
    const scalar T1 = T0 + small;
    return (func_->value(T1) - func_->value(T0))/small;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::saturationPressureModels::function::pSat
(
    const volScalarField::Internal& T
) const
{
    tmp<volScalarField::Internal> tpSat
    (
        volScalarField::Internal::New
        (
            IOobject::groupName("pSat", T.group()),
            T.mesh(),
            dimensionedScalar(dimPressure, 0)
        )
    );

    volScalarField::Internal& pSat = tpSat.ref();

    pSat.primitiveFieldRef() = func_->value(T.primitiveField());

    return tpSat;
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::saturationPressureModels::function::derivative
(
    const volScalarField::Internal& T
) const
{
    NotImplemented;
    tmp<volScalarField::Internal> tpSat
    (
        volScalarField::Internal::New
        (
            IOobject::groupName("pSat", T.group()),
            T.mesh(),
            dimensionedScalar(dimPressure/dimTemperature, 0)
        )
    );

    return tpSat;
}


// ************************************************************************* //
