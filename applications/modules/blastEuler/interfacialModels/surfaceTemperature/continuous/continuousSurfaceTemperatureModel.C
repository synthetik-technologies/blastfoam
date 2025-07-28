/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024-2025
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "continuousSurfaceTemperatureModel.H"
#include "phaseSystem.H"
#include "orderedPhasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceTemperatureModels
{
    defineTypeNameAndDebug(continuous, 0);
    addToRunTimeSelectionTable(surfaceTemperatureModel, continuous, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceTemperatureModels::continuous::continuous
(
    const dictionary& dict,
    const phaseModel& phase
)
:
    surfaceTemperatureModel(dict, phase)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceTemperatureModels::continuous::~continuous()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::surfaceTemperatureModels::continuous::Ts() const
{
    volScalarField alphaSum
    (
        volScalarField::New
        (
            "alphaSum",
            phase_.mesh(),
            0.0
        )
    );
    volScalarField alphaTSum
    (
        volScalarField::New
        (
            "alphaTSum",
            phase_.mesh(),
            dimensionedScalar(dimTemperature, 0.0)
        )
    );
    const phaseSystem& fluid = phase_.fluid();
    const volScalarField& T = phase_.T();
    forAllConstIter(phaseSystem::phasePairTable, fluid.phasePairs(), iter)
    {
        if
        (
            isA<orderedPhasePair>(iter()())
         && (&(iter()->dispersed())) == (&phase_)
        )
        {
            const phasePair& pair = iter()();
            const phaseModel& other = pair.continuous();

            alphaSum += other;
            alphaTSum += 0.5*other*(other.T() + T);
        }
    }


    return volScalarField::New
    (
        IOobject::groupName("Tsurface", phase_.name()),
        alphaTSum/max(alphaSum, phase_.residualAlpha())
    );
}


// ************************************************************************* //
