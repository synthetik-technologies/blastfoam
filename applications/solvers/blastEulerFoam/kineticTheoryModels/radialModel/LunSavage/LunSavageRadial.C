/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "LunSavageRadial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace radialModels
{
    defineTypeNameAndDebug(LunSavage, 0);

    addToRunTimeSelectionTable
    (
        radialModel,
        LunSavage,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::LunSavage::LunSavage
(
    const dictionary& dict,
    const kineticTheorySystem& kt
)
:
    radialModel(dict, kt)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::LunSavage::~LunSavage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::LunSavage::gs0
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    if (&phase1 != &phase2)
    {
        return
            volScalarField::New
            (
                "gs0prime",
                phase1.mesh(),
                dimensionedScalar(dimless, 0.0)
            );
    }
    return pow(1 - phase1/phase1.alphaMax(), -2.5*phase1.alphaMax());
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::LunSavage::gs0prime
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    if (&phase1 != &phase2)
    {
        return
            volScalarField::New
            (
                "gs0prime",
                phase1.mesh(),
                dimensionedScalar(dimless, 0.0)
            );
    }
    return 2.5*pow(1 - phase1/phase1.alphaMax(), -2.5*phase1.alphaMax() - 1);
}


// ************************************************************************* //
