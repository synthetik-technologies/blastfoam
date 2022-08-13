/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Synthetik Applied Technologies: |    Calculate mach number with
                                                blastFoam thermo
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

#include "blastMachNo.H"
#include "fluidBlastThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(blastMachNo, 0);
    addToRunTimeSelectionTable(functionObject, blastMachNo, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::blastMachNo::calc()
{
    if (!foundObject<volVectorField>(UName_))
    {
        return false;
    }
    const volVectorField& U = lookupObject<volVectorField>(UName_);

    if
    (
        foundObject<fluidBlastThermo>
        (
            IOobject::groupName(fluidBlastThermo::typeName, phaseName_)
        )
    )
    {
        tmp<volScalarField> speedOfSound
        (
            lookupObject<fluidBlastThermo>
            (
                IOobject::groupName(fluidBlastThermo::typeName, phaseName_)
            ).speedOfSound()
        );

        return store
        (
            resultName_,
            mag(U)/max(speedOfSound, dimensionedScalar(dimVelocity, small))
        );
    }
    else if
    (
        foundObject<fluidThermo>
        (
            IOobject::groupName(fluidThermo::typeName, phaseName_)
        )
    )
    {
        const fluidThermo& thermo
        (
            lookupObject<fluidThermo>
            (
                IOobject::groupName(fluidThermo::typeName, phaseName_)
            )
        );
        tmp<volScalarField> speedOfSound
        (
            sqrt(thermo.Cp()/thermo.Cv()/thermo.psi())
        );

        return store
        (
            resultName_,
            mag(U)/max(speedOfSound, dimensionedScalar(dimVelocity, small))
        );
    }
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::blastMachNo::blastMachNo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression
    (
        name,
        runTime,
        dict,
        IOobject::groupName("Ma", dict.lookupOrDefault("phaseName", word::null)),
        IOobject::groupName("U", dict.lookupOrDefault("phaseName", word::null))
    ),
    phaseName_(dict.lookupOrDefault("phaseName", word::null)),
    systemName_
    (
        dict.lookupOrDefault
        (
            "systemName",
            IOobject::groupName(basicThermo::dictName, phaseName_)
        )
    ),
    UName_(dict.lookupOrDefault("UName", IOobject::groupName("U", phaseName_)))
{
    if (!dict.lookupOrDefault("executeAtStart", false))
    {
        executeAtStart_ = false;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::blastMachNo::~blastMachNo()
{}


// ************************************************************************* //
