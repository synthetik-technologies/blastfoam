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

#ifndef blastCompressibleTurbulenceModel_H
#define blastCompressibleTurbulenceModel_H


#include "blastCompressibleTurbulenceModelFwd.H"
#include "CompressibleTurbulenceModel.H"
#include "ThermalDiffusivity.H"
#include "EddyDiffusivity.H"
#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"
#include "fluidThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef ThermalDiffusivity<CompressibleTurbulenceModel<fluidThermoModel>>
        blastCompressibleTurbulenceModel;

    namespace blast
    {
        typedef laminarModel<blastCompressibleTurbulenceModel> laminarModel;
        typedef RASModel<EddyDiffusivity<blastCompressibleTurbulenceModel>> RASModel;
        typedef LESModel<EddyDiffusivity<blastCompressibleTurbulenceModel>> LESModel;
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
