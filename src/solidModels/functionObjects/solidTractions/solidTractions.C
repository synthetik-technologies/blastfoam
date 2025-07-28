/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*----------------------------------------------------------------------------*/

#include "solidTractions.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(solidTractions, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        solidTractions,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::solidTractions::solidTractions
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, t, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::solidTractions::read(const dictionary& dict)
{
    return fvMeshFunctionObject::read(dict);
}


bool Foam::functionObjects::solidTractions::execute()
{
    // Lookup the stress field
    tmp<surfaceSymmTensorField> tsigmaf;
    if (mesh_.foundObject<surfaceSymmTensorField>("sigmaf"))
    {
        tsigmaf = tmp<surfaceSymmTensorField>
        (
            mesh_.lookupObject<surfaceSymmTensorField>("sigmaf")
        );
    }
    else
    {
        tsigmaf = fvc::interpolate
        (
            mesh_.lookupObject<volSymmTensorField>("sigma")
        );
    }
    const surfaceSymmTensorField& sigmaf = tsigmaf();

    if (mesh_.foundObject<volTensorField>("Ff"))
    {
        // Total Lagrangian (face based)
        // The mesh is in its initial configuration

        // Lookup the inverse deformation gradient
        const surfaceTensorField& Ffinv =
            mesh_.lookupObject<surfaceTensorField>("Ffinv");

        return store
        (
            "traction",
            (
                (Ffinv.T() & (mesh_.Sf()/mesh_.magSf()))
              & sigmaf
            )
        );
    }
    else if (mesh_.foundObject<volTensorField>("F"))
    {
        // Total Lagrangian
        // The mesh is in its initial configuration

        // Lookup the inverse deformation gradient
        const volTensorField& Finv =
            mesh_.lookupObject<volTensorField>("Finv");

        return store
        (
            "traction",
            (
                (fvc::interpolate(Finv.T()) & (mesh_.Sf()/mesh_.magSf()))
              & sigmaf
            )
        );
    }
    else
    {
        return store
        (
            "traction",
            ((mesh_.Sf()/mesh_.magSf()) & sigmaf)
        );
    }
}



bool Foam::functionObjects::solidTractions::write()
{
    return writeObject("traction");
}

// ************************************************************************* //
