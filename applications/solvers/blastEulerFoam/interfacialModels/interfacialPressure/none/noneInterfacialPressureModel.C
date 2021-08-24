/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "noneInterfacialPressureModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interfacialPressureModels
{
    defineTypeNameAndDebug(none, 0);
    addToRunTimeSelectionTable(interfacialPressureModel, none, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfacialPressureModels::none::none
(
    const dictionary& dict,
    const phasePair& pair
)
:
    interfacialPressureModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfacialPressureModels::none::~none()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::interfacialPressureModels::none::PI() const
{
    return volScalarField::New
    (
        "none:PI",
        this->pair_.phase1().mesh(),
        dimensionedScalar(dimPressure, 0.0)
    );
}


Foam::scalar
Foam::interfacialPressureModels::none::cellPI(const label celli) const
{
    return 0.0;
}


Foam::scalar
Foam::interfacialPressureModels::none::celldPIdAlpha
(
    const label celli,
    const label phasei
) const
{
    return 0.0;
}


Foam::scalar
Foam::interfacialPressureModels::none::celldPIde
(
    const label celli,
    const label phasei
) const
{
    return 0.0;
}


// ************************************************************************* //
