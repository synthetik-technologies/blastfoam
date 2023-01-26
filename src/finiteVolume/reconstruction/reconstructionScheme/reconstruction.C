/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020
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
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //
Foam::word Foam::reconstruction::scheme(const word& name)
{
    return "reconstruct(" + name + ")";
}


Foam::word Foam::reconstruction::scheme
(
    const word& name,
    const fvMesh& mesh,
    const bool fail,
    const bool overwrite
)
{
    return scheme(name, word::null, mesh, fail, overwrite);
}


Foam::word Foam::reconstruction::scheme
(
    const word& baseName,
    const word& phaseName,
    const fvMesh& mesh,
    const bool fail,
    const bool overwrite
)
{
    const word name(IOobject::groupName(baseName, phaseName));
    word baseScheme(scheme(baseName));
    word nameScheme(scheme(name));

    if (mesh.schemesDict().subDict("interpolationSchemes").found(nameScheme))
    {
        return nameScheme;
    }
    else if (mesh.schemesDict().subDict("interpolationSchemes").found(baseScheme))
    {
        return baseScheme;
    }
    else if (fail && overwrite)
    {
        FatalErrorInFunction
            << "Riemann fluxes are used, but no limiter is " << nl
            << "specified for " << name << "." << nl
            << "Please specify " << string(nameScheme)
            << " or " << string(baseScheme) << endl
            << "This may result in unstable solutions." << endl
            << abort(FatalError);
    }
    else if (overwrite)
    {
        WarningInFunction
            << "Riemann fluxes are used, but no limiter is " << nl
            << "specified for " << name << "." << nl
            << "This may result in unstable solutions." << nl
            << "Please specify " << string(nameScheme)
            << " or " << string(baseScheme) << endl;

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
