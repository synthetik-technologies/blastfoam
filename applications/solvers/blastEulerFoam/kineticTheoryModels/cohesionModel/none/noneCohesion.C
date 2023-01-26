/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "noneCohesion.H"
#include "kineticTheoryModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace kineticTheoryModels
{
namespace cohesionModels
{
    defineTypeNameAndDebug(none, 0);

    addToRunTimeSelectionTable
    (
        cohesionModel,
        none,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::cohesionModels::none::none
(
    const dictionary& dict,
    const kineticTheoryModel& kt
)
:
    cohesionModel(dict, kt)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModels::cohesionModels::none::~none()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::cohesionModels::none::nu() const
{
    return
        volScalarField::New
        (
            "none:PsCoh",
            kt_.phase().mesh(),
            dimensionedScalar(dimViscosity, 0)
        );
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::cohesionModels::none::Ps() const
{
    return
        volScalarField::New
        (
            "none:PsCoh",
            kt_.phase().mesh(),
            dimensionedScalar(dimPressure, 0)
        );
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::cohesionModels::none::dPsdAlpha() const
{
    return volScalarField::New
        (
            "none:dPsdAlpha",
            kt_.phase().mesh(),
            dimensionedScalar(dimPressure, 0)
        );
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::cohesionModels::none::dPsdTheta() const
{
    return
        volScalarField::New
        (
            "none:dPsdTheta",
            kt_.phase().mesh(),
            dimensionedScalar(dimPressure/sqr(dimVelocity), 0.0)
        );
}


Foam::tmp<Foam::volScalarField>
Foam::kineticTheoryModels::cohesionModels::none::dissipationSource
(
    const dimensionedScalar&
) const
{
    return
        volScalarField::New
        (
            "none:dissipationSource",
            kt_.phase().mesh(),
            dimensionedScalar(dimDensity*sqr(dimVelocity), 0.0)
        );
}


// ************************************************************************* //
