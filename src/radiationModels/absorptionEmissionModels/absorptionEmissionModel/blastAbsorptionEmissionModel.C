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

#include "blastAbsorptionEmissionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiationModels
    {
        defineTypeNameAndDebug(blastAbsorptionEmissionModel, 0);
        defineRunTimeSelectionTable(blastAbsorptionEmissionModel, dictionary);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::blastAbsorptionEmissionModel::blastAbsorptionEmissionModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh)
{}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiationModels::blastAbsorptionEmissionModel::~blastAbsorptionEmissionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::radiationModels::blastAbsorptionEmissionModel::ai
(
    const label celli,
    const label bandI
) const
{
    return aDispi(celli, bandI) + aConti(celli, bandI);
}


Foam::scalar
Foam::radiationModels::blastAbsorptionEmissionModel::aConti
(
    const label celli,
    const label bandI
) const
{
    return 0.0;
}


Foam::scalar
Foam::radiationModels::blastAbsorptionEmissionModel::aDispi
(
    const label celli,
    const label bandI
) const
{
    return 0.0;
}


Foam::scalar
Foam::radiationModels::blastAbsorptionEmissionModel::ei
(
    const label celli,
    const label bandI
) const
{
    return eDispi(celli, bandI) + eConti(celli, bandI);
}


Foam::scalar
Foam::radiationModels::blastAbsorptionEmissionModel::eConti
(
    const label celli,
    const label bandI
) const
{
    return 0.0;
}


Foam::scalar
Foam::radiationModels::blastAbsorptionEmissionModel::eDispi
(
    const label celli,
    const label bandI
) const
{
    return 0.0;
}


Foam::scalar
Foam::radiationModels::blastAbsorptionEmissionModel::Ei
(
    const label celli,
    const label bandI
) const
{
    return EDispi(celli, bandI) + EConti(celli, bandI);
}


Foam::scalar
Foam::radiationModels::blastAbsorptionEmissionModel::EConti
(
    const label celli,
    const label bandI
) const
{
    return 0.0;
}


Foam::scalar
Foam::radiationModels::blastAbsorptionEmissionModel::EDispi
(
    const label celli,
    const label bandI
) const
{
    return 0.0;
}
//

// ************************************************************************* //
