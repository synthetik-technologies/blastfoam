/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2024
     \\/     M anipulation  | Synthetik Applied Technology
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

#include "reconstruction.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(reconstruction, 0);

    wordHashSet* reconstruction::limiterFreeSchemesPtr = nullptr;

    wordHashSet& reconstruction::limiterFreeSchemes()
    {
        if (!limiterFreeSchemesPtr)
        {
            limiterFreeSchemesPtr = new wordHashSet();
        }
        return *limiterFreeSchemesPtr;
    }

    void reconstruction::destroyLimiterFreeSchemes()
    {
        if (limiterFreeSchemesPtr)
        {
            deleteDemandDrivenData(limiterFreeSchemesPtr);
        }
    }
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //
Foam::word Foam::reconstruction::scheme(const word& name)
{
    return "reconstruct(" + name + ")";
}


Foam::word Foam::reconstruction::scheme
(
    const word& name,
    const word& type,
    const fvSchemes& schemes,
    const bool fail,
    const bool overwrite
)
{
    return scheme
    (
        IOobject::member(name),
        IOobject::group(name),
        type,
        schemes,
        fail,
        overwrite
    );
}


Foam::word Foam::reconstruction::scheme
(
    const word& baseName,
    const word& phaseName,
    const word& type,
    const fvSchemes& schemes,
    const bool fail,
    const bool overwrite
)
{
    const dictionary& interpDict =
        schemes.dict().subDict("interpolationSchemes");
    const word name(IOobject::groupName(baseName, phaseName));
    word baseScheme(scheme(baseName));
    word nameScheme(scheme(name));
    word typeScheme(scheme(type));

    // Exact match, no pattern
    if (interpDict.found(nameScheme, false, false))
    {
        return nameScheme;
    }
    else if (interpDict.found(baseScheme, false, false))
    {
        return baseScheme;
    }
    else if (interpDict.found(typeScheme, false, false))
    {
        return typeScheme;
    }

    // Patterns allowed
    else if (interpDict.found(nameScheme))
    {
        return nameScheme;
    }
    else if (interpDict.found(baseScheme))
    {
        return baseScheme;
    }
    else if (interpDict.found(typeScheme))
    {
        return typeScheme;
    }

    // Default
    else if (interpDict.found("defaultReconstruction"))
    {
        return "defaultReconstruction";
    }

    // Not found
    else if (fail && overwrite)
    {
        FatalErrorInFunction
            << "Riemann fluxes are used, but no limiter is " << nl
            << "specified for " << name << "." << nl
            << "Please specify " << string(baseScheme)
            << ", " << string(nameScheme)
            << ", or " << string(typeScheme) << endl
            << "This may result in unstable solutions." << endl
            << abort(FatalError);
    }
    else if (overwrite)
    {
        WarningInFunction
            << "Riemann fluxes are used, but no limiter is " << nl
            << "specified for " << name << "." << nl
            << "This may result in unstable solutions." << nl
            << "Please specify " << string(baseScheme)
            << ", " << string(nameScheme)
            << ", or " << string(typeScheme) << endl;

    }
    return nameScheme;
}

Foam::word Foam::reconstruction::ownName(const word& name)
{
    return
        IOobject::groupName
        (
            IOobject::member(name) + "Own",
            IOobject::group(name)
        );
}


Foam::word Foam::reconstruction::neiName(const word& name)
{
    return
        IOobject::groupName
        (
            IOobject::member(name) + "Nei",
            IOobject::group(name)
        );
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::reconstruction::reconstruction()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::reconstruction::~reconstruction()
{}


// ************************************************************************* //
