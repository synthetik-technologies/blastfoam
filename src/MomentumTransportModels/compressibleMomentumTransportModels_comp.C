/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2025-06-09 Jeff Heylmun     : Compile compresibility correction momentum
                              transport models
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

#include "makeCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

#include "kOmega_comp.H"
makeRASModel(kOmega_comp);

#include "kOmega2006_comp.H"
makeRASModel(kOmega2006_comp);

#include "kOmegaSST_comp.H"
makeRASModel(kOmegaSST_comp);

#include "kOmegaSSTSAS_comp.H"
makeRASModel(kOmegaSSTSAS_comp);

#include "kOmegaSSTLM_comp.H"
makeRASModel(kOmegaSSTLM_comp);


// -------------------------------------------------------------------------- //
// LES models
// -------------------------------------------------------------------------- //

#include "kOmegaSSTDES_comp.H"
makeLESModel(kOmegaSSTDES_comp);


// ************************************************************************* //
