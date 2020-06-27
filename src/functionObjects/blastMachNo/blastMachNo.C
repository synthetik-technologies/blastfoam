/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Jeff Heylmun:    | Calculate mach number with blastFoam thermo
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CompressibleSystem>
Foam::functionObjects::blastMachNo<CompressibleSystem>::blastMachNo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    phaseName_(dict.lookupOrDefault("phaseName", word::null)),
    systemName_(IOobject::groupName(systemType_, phaseName_)),
    resultName_(IOobject::groupName("Ma", phaseName_)),
    UName_(IOobject::groupName("U", phaseName_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CompressibleSystem>
Foam::functionObjects::blastMachNo<CompressibleSystem>::~blastMachNo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompressibleSystem>
bool Foam::functionObjects::blastMachNo<CompressibleSystem>::read
(
    const dictionary& dict
)
{
    return true;
}


template<class CompressibleSystem>
bool Foam::functionObjects::blastMachNo<CompressibleSystem>::execute()
{
    if
    (
        foundObject<volVectorField>(UName_)
     && foundObject<CompressibleSystem>(systemName_)
    )
    {
        tmp<volScalarField> speedOfSound
        (
            lookupObject<CompressibleSystem>(systemName_).speedOfSound()
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


template<class CompressibleSystem>
bool Foam::functionObjects::blastMachNo<CompressibleSystem>::write()
{
    writeObject(resultName_);
    return true;
}


template<class CompressibleSystem>
bool Foam::functionObjects::blastMachNo<CompressibleSystem>::clear()
{
    return clearObject(resultName_);
}

// ************************************************************************* //
