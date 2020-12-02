/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

InClass
    Foam::psiChemistryModel

Description
    Creates chemistry model instances templated on the type of thermodynamics

\*---------------------------------------------------------------------------*/

#include "makeChemistryModel.H"

#include "multicomponentFluidThermoTypes.H"

#include "StandardChemistryModel.H"
#include "TDACChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Make base types
    makeChemistryModelFromTypes
    (
        constTransport,
        eConst,
        idealGas
    );
    makeChemistryModelFromTypes
    (
        constTransport,
        hConst,
        idealGas
    );
    makeChemistryModelFromTypes
    (
        constTransport,
        janaf,
        idealGas
    );

    makeChemistryModelFromTypes
    (
        sutherlandTransport,
        eConst,
        idealGas
    );
    makeChemistryModelFromTypes
    (
        sutherlandTransport,
        hConst,
        idealGas
    );
    makeChemistryModelFromTypes
    (
        sutherlandTransport,
        janaf,
        idealGas
    );

    makeChemistryModelFromTypes
    (
        constTransport,
        eConst,
        perfectGas
    );
    makeChemistryModelFromTypes
    (
        constTransport,
        hConst,
        perfectGas
    );
    makeChemistryModelFromTypes
    (
        constTransport,
        janaf,
        perfectGas
    );

    makeChemistryModelFromTypes
    (
        sutherlandTransport,
        eConst,
        perfectGas
    );
    makeChemistryModelFromTypes
    (
        sutherlandTransport,
        hConst,
        perfectGas
    );
    makeChemistryModelFromTypes
    (
        sutherlandTransport,
        janaf,
        perfectGas
    );
}

// ************************************************************************* //
