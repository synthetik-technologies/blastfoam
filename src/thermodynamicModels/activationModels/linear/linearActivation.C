/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "linearActivation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace activationModels
{
    defineTypeNameAndDebug(linearActivation, 0);
    addToRunTimeSelectionTable(activationModel, linearActivation, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::activationModels::linearActivation::linearActivation
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    activationModel(mesh, dict, phaseName, true),
    vDet_("vDet", dimVelocity, dict),
    tIgn_
    (
        IOobject
        (
            IOobject::groupName("tIgn", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimTime, great)
    )
{
    // Are delays used to calculate ignition time
    // if no, the closest point is used
    Switch useDelay(dict.lookupOrDefault("delayOffset", false));
    scalarField tIgn(mesh.nCells(), great);
    forAll(this->detonationPoints_, pointi)
    {
        const detonationPoint& dp = this->detonationPoints_[pointi];
        forAll(tIgn_, celli)
        {
            scalar t =
                (useDelay ? dp.delay() : 0.0)
              + mag(this->mesh().C()[celli] - dp)/vDet_.value();
            if (t < tIgn[celli])
            {
                tIgn_[celli] = dp.delay() + t;
                tIgn[celli] = t;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::activationModels::linearActivation::~linearActivation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::activationModels::linearActivation::delta() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("linear:R", lambda_.group()),
                lambda_.time().timeName(),
                lambda_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            lambda_.mesh(),
            dimensionedScalar
            (
                "delta",
                inv(dimTime),
                0.0
            )
        )
    );
}


void Foam::activationModels::linearActivation::correct()
{
    if (min(lambda_.oldTime()).value() == 1)
    {
        return;
    }
    volScalarField::Internal lambda0
    (
        pos0(lambda_.time() - lambda_.time().deltaT() - tIgn_)
    );
    volScalarField::Internal diff(pos0(lambda_.time() - tIgn_) - lambda0);
    lambda_.ref() = max(lambda0 + diff*(this->f() - this->f0()), lambda_());
}

// ************************************************************************* //
