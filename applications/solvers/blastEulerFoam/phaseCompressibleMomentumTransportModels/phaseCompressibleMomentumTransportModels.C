/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2019 OpenFOAM Foundation
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

#include "phaseDynamicMomentumTransportModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Stokes.H"
makeLaminarModel(Stokes);

#include "generalisedNewtonian.H"
makeLaminarModel(generalisedNewtonian);

#include "kEpsilon.H"
makeRASModel(kEpsilon);

// #include "LaheyKEpsilon.H"
// makeRASModel(LaheyKEpsilon);

// #include "kOmegaSST.H"
// makeRASModel(kOmegaSST);

// #include "kOmegaSSTSato.H"
// makeRASModel(kOmegaSSTSato);

// #include "continuousGasKEpsilon.H"
// makeRASModel(continuousGasKEpsilon);

// #include "mixtureKEpsilon.H"
// makeRASModel(mixtureKEpsilon);

#include "Smagorinsky.H"
makeLESModel(Smagorinsky);

#include "kEqn.H"
makeLESModel(kEqn);

// #include "SmagorinskyZhang.H"
// makeLESModel(SmagorinskyZhang);

// #include "NicenoKEqn.H"
// makeLESModel(NicenoKEqn);

// #include "continuousGasKEqn.H"
// makeLESModel(continuousGasKEqn);


// ************************************************************************* //
