/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
Typedef
    Foam::fiveEqnCompressibleTurbulenceModel
Description
    Typedef for phaseCompressibleTurbulenceModel
\*---------------------------------------------------------------------------*/

#include "CompressibleMomentumTransportModel.H"
#include "fluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"
#include "makeBlastMomentumTransportModel.H"

#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeMomentumTransportModelTypes
(
    geometricOneField,
    volScalarField,
    compressibleMomentumTransportModel,
    CompressibleMomentumTransportModel,
    fluidThermoModel
);


#define makeLaminarModel(Type)                                                 \
    makeTemplatedLaminarModel (fluidThermoModelCompressibleMomentumTransportModel, laminar, Type)

#define makeRASModel(Type)                                                     \
    makeTemplatedMomentumTransportModel                                        \
    (fluidThermoModelCompressibleMomentumTransportModel, RAS, Type)

#define makeLESModel(Type)                                                     \
    makeTemplatedMomentumTransportModel                                        \
    (fluidThermoModelCompressibleMomentumTransportModel, LES, Type)


// ************************************************************************* //
