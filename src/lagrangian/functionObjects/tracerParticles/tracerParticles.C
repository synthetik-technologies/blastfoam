/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "tracerParticles.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(tracerParticles, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        tracerParticles,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::tracerParticles::tracerParticles
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    g_
    (
        IOobject
        (
            "g",
            time_.constant(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedVector(dimAcceleration, Zero)
    ),
    rho_
    (
        mesh_.lookupObject<volScalarField>
        (
            dict.lookupOrDefault<word>("rho", "rho")
        )
    ),
    U_
    (
        mesh_.lookupObject<volVectorField>
        (
            dict.lookupOrDefault<word>("U", "U")
        )
    ),
    mu_
    (
        mesh_.lookupObject<volScalarField>
        (
            IOobject::groupName("thermo:mu", U_.group())
        )
    ),
    kinematicCloudName_
    (
        dict.lookupOrDefault<word>("kinematicCloud", "kinematicCloud")
    ),
    kinematicCloud_
    (
        kinematicCloudName_,
        rho_,
        U_,
        mu_,
        g_
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::tracerParticles::~tracerParticles()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::tracerParticles::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);

    return true;
}


bool Foam::functionObjects::tracerParticles::execute()
{
    kinematicCloud_.evolve();

    return true;
}


bool Foam::functionObjects::tracerParticles::write()
{
    return true;
}


// ************************************************************************* //
