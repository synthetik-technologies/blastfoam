/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "intact.H"
#include "burstPolyPatchBase.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(intact, 0);
    addToRunTimeSelectionTable(functionObject, intact, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::intact::intact
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::intact::~intact()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::functionObjects::intact::write()
{
    if (this->mesh_.time().timeIndex() > 0)
    {
        volScalarField intact
        (
            IOobject
            (
                "intact",
                obr_.time().timeName(),
                obr_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh_,
            dimensionedScalar("0", dimless, 0.0)
        );

        forAll(this->mesh_.boundaryMesh(), patchi)
        {
            const polyPatch& pp = this->mesh_.boundaryMesh()[patchi];
            if (isA<burstPolyPatchBase>(pp))
            {
                const burstPolyPatchBase& bppb =
                    dynamicCast<const burstPolyPatchBase>(pp);
                intact.boundaryFieldRef()[patchi] == bppb.intact();
            }
            else if (isA<wallPolyPatch>(pp))
            {
                intact.boundaryFieldRef() = 1.0;
            }
        }
        return intact.write();
    }
    return false;
}


// ************************************************************************* //
