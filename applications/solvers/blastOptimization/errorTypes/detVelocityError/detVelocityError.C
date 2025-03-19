/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technologies
     \\/     M anipulation  |
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

#include "detVelocityError.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace errorTypes
{
    defineTypeNameAndDebug(detVelocity, 0);
    addToRunTimeSelectionTable(errorType, detVelocity, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::errorTypes::detVelocity::detVelocity
(
    const Time& runTime,
    const dictionary& dict,
    const word& region
)
:
    errorType(runTime, dict, region, LATEST),
    pName_(dict.lookupOrDefault<word>("pName", "p")),
    point1_(readConfigValue<vector>("point1", dict)),
    point2_(readConfigValue<vector>("point2", dict)),
    startTime_(-1.0),
    endTime_(-1.0),
    stopWhenReached_(dict.lookupOrDefault<Switch>("stopWhenReached", true)),
    pTarget_(readConfigValue<scalar>("pTarget", dict))
{
    value_ = 0.0;
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::errorTypes::detVelocity::~detVelocity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::errorTypes::detVelocity::update()
{
    const fvMesh& mesh = this->mesh();
    const Time& runTime = mesh.time();

    // Look up pressure
    const volScalarField& p
    (
        mesh.lookupObject<volScalarField>(pName_)
    );

    if (startTime_ < 0)
    {
        label celli = mesh.findCell(point1_);
        if (celli >= 0 && p[celli] > pTarget_)
        {
            point1_ = mesh.C()[celli];
            startTime_ = runTime.value();
        }
        reduce(startTime_, maxOp<scalar>());

        if (Pstream::parRun() && startTime_ > 0)
        {
            if (celli < 0)
            {
                point1_ = -great*vector::one;
            }
            reduce(point1_, maxOp<vector>());
        }
    }
    if (endTime_ < 0)
    {
        label celli = mesh.findCell(point2_);
        if (celli >= 0 && p[celli] > pTarget_)
        {
            point2_ = mesh.C()[celli];
            endTime_ = runTime.value();
        }
        reduce(endTime_, maxOp<scalar>());

        if (Pstream::parRun() && endTime_ > 0)
        {
            if (celli < 0)
            {
                point2_ = -great*vector::one;
            }
            reduce(point2_, maxOp<vector>());
        }
    }

    if (endTime_ > 0 && startTime_ > 0)
    {
        if (stopWhenReached_ && value_ != 0)
        {
            const_cast<Time&>(runTime).stopAt
            (
                Time::stopAtControl::nextWrite
            );
        }
        value_ = mag(point2_ - point1_)/(endTime_ - startTime_);
    }
}

// ************************************************************************* //
