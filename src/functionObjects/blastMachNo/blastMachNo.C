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
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::blastMachNo::blastMachNo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    phaseName_(dict.lookupOrDefault("phaseName", word::null)),
    systemName_(IOobject::groupName(basicThermo::dictName, phaseName_)),
    resultName_(IOobject::groupName("Ma", phaseName_)),
    UName_(IOobject::groupName("U", phaseName_))
{
    if (!dict.lookupOrDefault("executeAtStart", false))
    {
        executeAtStart_ = false;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::blastMachNo::~blastMachNo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::blastMachNo::read
(
    const dictionary& dict
)
{
    return true;
}


bool Foam::functionObjects::blastMachNo::execute()
{
    if
    (
        foundObject<volVectorField>(UName_)
     && foundObject<fluidBlastThermo>(systemName_)
    )
    {
        tmp<volScalarField> speedOfSound
        (
            lookupObject<fluidBlastThermo>(systemName_).speedOfSound()
        );
        speedOfSound.ref().max(small);

        const volVectorField& U = lookupObject<volVectorField>(UName_);

        return store
        (
            resultName_,
            mag(U)/speedOfSound
        );
    }
    else
    {
        Warning
            << "    functionObjects::" << type() << " " << name()
            << " failed to execute." << endl;

        return false;
    }
}


bool Foam::functionObjects::blastMachNo::write()
{
    writeObject(resultName_);
    return true;
}


bool Foam::functionObjects::blastMachNo::clear()
{
    return clearObject(resultName_);
}

// ************************************************************************* //
