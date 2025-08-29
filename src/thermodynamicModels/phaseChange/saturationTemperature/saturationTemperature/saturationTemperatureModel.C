/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
28-08-2025  Synthetik Applied Technologies: Added single input functions
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

#include "saturationTemperatureModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(saturationTemperatureModel, 0);
    defineRunTimeSelectionTable(saturationTemperatureModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationTemperatureModel::saturationTemperatureModel()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationTemperatureModel::~saturationTemperatureModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::saturationTemperatureModel::Tsat
(
    const volScalarField::Internal& p
) const
{
    tmp<volScalarField::Internal> tTsat
    (
        volScalarField::Internal::New
        (
            "Tsat",
            p.mesh(),
            dimensionedScalar(dimTemperature, 0)
        )
    );

    volScalarField::Internal& Tsat = tTsat.ref();

    forAll(Tsat, celli)
    {
        Tsat[celli] = this->Tsat(p[celli]);
    }

    return tTsat;
}


Foam::tmp<Foam::volScalarField> Foam::saturationTemperatureModel::Tsat
(
    const volScalarField& p
) const
{
    tmp<volScalarField> tTsat
    (
        volScalarField::New
        (
            "Tsat",
            p.mesh(),
            dimensionedScalar(dimTemperature, 0)
        )
    );

    volScalarField& Tsat = tTsat.ref();

    forAll(Tsat, celli)
    {
        Tsat[celli] = this->Tsat(p[celli]);
    }

    volScalarField::Boundary& TsatBf = Tsat.boundaryFieldRef();

    forAll(Tsat.boundaryField(), patchi)
    {
        scalarField& Tsatp = TsatBf[patchi];
        const scalarField& pp = p.boundaryField()[patchi];

        forAll(Tsatp, facei)
        {
            Tsatp[facei] = this->Tsat(pp[facei]);
        }
    }

    return tTsat;
}

// ************************************************************************* //
