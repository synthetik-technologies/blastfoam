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

#include "LunPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace granularPressureModels
{
    defineTypeNameAndDebug(Lun, 0);

    addToRunTimeSelectionTable
    (
        granularPressureModel,
        Lun,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::granularPressureModels::Lun::Lun
(
    const dictionary& dict,
    const kineticTheorySystem& kt
)
:
    granularPressureModel(dict, kt)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::granularPressureModels::Lun::~Lun()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::Lun::granularPressure
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const volScalarField& Theta1,
    const volScalarField& Theta2,
    const volScalarField& g0,
    const dimensionedScalar& e
) const
{
    dimensionedScalar eta = 0.5*(1.0 + e);
    return phase1*phase1.rho()*Theta1*4.0*eta*phase2*g0;
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::Lun::
granularPressureByAlpha
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const volScalarField& Theta1,
    const volScalarField& Theta2,
    const volScalarField& g0,
    const volScalarField& g0Prime,
    const dimensionedScalar& e
) const
{
    dimensionedScalar eta = 0.5*(1.0 + e);

    if (&phase1 == &phase2)
    {
        return 4.0*phase1*phase1.rho()*Theta1*eta*(2.0*g0 + phase1*g0Prime);
    }
    return 4.0*phase2*phase1.rho()*Theta1*eta*(g0 + phase1*g0Prime);
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::granularPressureModels::Lun::
granularPressureByTheta
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const volScalarField& Theta1,
    const volScalarField& Theta2,
    const volScalarField& g0,
    const dimensionedScalar& e
) const
{
    dimensionedScalar eta = 0.5*(1.0 + e);
    return phase1*phase1.rho()*4.0*eta*phase2*g0;
}


// ************************************************************************* //
