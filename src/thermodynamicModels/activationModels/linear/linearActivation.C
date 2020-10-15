/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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
    activationModel(mesh, dict, phaseName),
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
    scalarField distance(tIgn_.size(), great);
    forAll(this->detonationPoints_, pointi)
    {
        const detonationPoint& dp = this->detonationPoints_[pointi];
        forAll(distance, celli)
        {
            scalar d = mag(mesh.cellCentres()[celli] - dp);
            if (d < distance[celli])
            {
                distance[celli] = d;
                tIgn_[celli] =
                    dp.delay() + mag(mesh_.C()[celli] - dp)/vDet_.value();
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
    tmp<volScalarField> deltaLambdaTmp
    (
        new volScalarField
        (
            IOobject
            (
                "deltaLambda",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", inv(dimTime), 0.0)
        )
    );
    volScalarField& deltaLambda = deltaLambdaTmp.ref();

    dimensionedScalar dt(lambda_.time().deltaT());
    dimensionedScalar t(lambda_.time());

    forAll(lambda_, celli)
    {
        if (t.value() > tIgn_[celli])
        {
            deltaLambda[celli] = (1.0 - lambda_.oldTime()[celli])/dt.value();
        }
    }

    return deltaLambdaTmp;
}

// ************************************************************************* //
