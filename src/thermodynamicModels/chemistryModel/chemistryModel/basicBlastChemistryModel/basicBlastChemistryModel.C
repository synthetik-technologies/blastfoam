/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "basicBlastChemistryModel.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(basicBlastChemistryModel, 0);
    defineRunTimeSelectionTable(basicBlastChemistryModel, thermo);
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::basicBlastChemistryModel::correct()
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicBlastChemistryModel::basicBlastChemistryModel
(
    const blastThermo& mixture
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName("chemistryProperties", mixture.phaseName()),
            mixture.mesh().time().constant(),
            mixture.mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mixture.mesh()),
    mixture_(dynamicCast<const multicomponentBlastThermo>(mixture)),
    thermo_(refCast<const blastThermo>(mixture)),
    chemistry_(lookup("chemistry")),
    deltaTChemIni_(lookup<scalar>("initialChemicalTimeStep")),
    deltaTChemMax_(lookupOrDefault("maxChemicalTimeStep", great)),
    deltaTChem_
    (
        IOobject
        (
            thermo_.phasePropertyName("deltaTChem"),
            mesh().time().constant(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedScalar(dimTime, deltaTChemIni_)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicBlastChemistryModel::~basicBlastChemistryModel()
{}


// ************************************************************************* //
