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

#include "CarnahanStarlingRadial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace radialModels
{
    defineTypeNameAndDebug(CarnahanStarling, 0);

    addToRunTimeSelectionTable
    (
        radialModel,
        CarnahanStarling,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::CarnahanStarling::CarnahanStarling
(
    const dictionary& dict,
    const masterSystem& system
)
:
    radialModel(dict, system)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::radialModels::CarnahanStarling::~CarnahanStarling()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::CarnahanStarling::gs0
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
                "gs0",
                phase1.mesh(),
                dimensionedScalar(dimless, 0.0)
            );
    }
    const volScalarField aBounded(min(phase1, phase1.alphaMax() - 1e-4));
    const volScalarField oneMinusAlpha(max(1.0 - aBounded, dimensionedScalar(dimless, 1e-4)));
    return
        1.0/oneMinusAlpha
      + 3.0*aBounded/(2.0*sqr(oneMinusAlpha))
      + sqr(aBounded)/(2.0*pow3(oneMinusAlpha));
}


Foam::scalar
Foam::kineticTheoryModels::radialModels::CarnahanStarling::cellgs0
(
    const label celli,
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    if (&phase1 != &phase2)
    {
        return 0.0;
    }
    const scalar aBounded = min(phase1[celli], phase1.alphaMax() - 1e-4);
    const scalar oneMinusAlpha = max(1.0 - aBounded, 1e-4);
    return
        1.0/oneMinusAlpha
      + 3.0*aBounded/(2.0*sqr(oneMinusAlpha))
      + sqr(aBounded)/(2.0*pow3(oneMinusAlpha));
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::radialModels::CarnahanStarling::gs0prime
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
    const volScalarField aBounded(min(phase1, phase1.alphaMax() - 1e-4));
    const volScalarField oneMinusAlpha(max(1.0 - aBounded, dimensionedScalar(dimless, 1e-4)));
    return
        2.5/sqr(oneMinusAlpha)
      + 4.0*aBounded/pow3(oneMinusAlpha)
      + 1.5*sqr(aBounded)/pow4(oneMinusAlpha);
}


Foam::scalar
Foam::kineticTheoryModels::radialModels::CarnahanStarling::cellgs0prime
(
    const label celli,
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    if (&phase1 != &phase2)
    {
        return 0.0;
    }
    const scalar aBounded = min(phase1[celli], phase1.alphaMax() - 1e-4);
    const scalar oneMinusAlpha = max(1.0 - aBounded, 1e-4);
    return
        2.5/sqr(oneMinusAlpha)
      + 4.0*aBounded/pow3(oneMinusAlpha)
      + 1.5*sqr(aBounded)/pow4(oneMinusAlpha);
}


// ************************************************************************* //
