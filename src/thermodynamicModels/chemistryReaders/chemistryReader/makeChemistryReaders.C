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

\*---------------------------------------------------------------------------*/

#include "fluidThermoModelTypes.H"

#include "chemistryReader.H"
#include "foamChemistryReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeChemistryReaderFromTypes(constTransport, eConst, idealGas);
makeChemistryReaderFromTypes(constTransport, hConst, idealGas);
makeChemistryReaderFromTypes(constTransport, janaf, idealGas);
makeChemistryReaderFromTypes(sutherlandTransport, eConst, idealGas);
makeChemistryReaderFromTypes(sutherlandTransport, hConst, idealGas);
makeChemistryReaderFromTypes(sutherlandTransport, janaf, idealGas);

makeChemistryReaderFromTypes(constTransport, eConst, perfectGas);
makeChemistryReaderFromTypes(constTransport, hConst, perfectGas);
makeChemistryReaderFromTypes(constTransport, janaf, perfectGas);
makeChemistryReaderFromTypes(sutherlandTransport, eConst, perfectGas);
makeChemistryReaderFromTypes(sutherlandTransport, hConst, perfectGas);
makeChemistryReaderFromTypes(sutherlandTransport, janaf, perfectGas);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
