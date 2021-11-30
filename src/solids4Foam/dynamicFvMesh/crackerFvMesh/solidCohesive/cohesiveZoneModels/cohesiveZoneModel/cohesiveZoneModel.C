/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    cohesiveZoneModel

\*---------------------------------------------------------------------------*/

#include "cohesiveZoneModel.H"
#include "volFields.H"
#include "fvc.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cohesiveZoneModel, 0);
    defineRunTimeSelectionTable(cohesiveZoneModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cohesiveZoneModel::cohesiveZoneModel
(
    const word& name,
    const fvPatch& patch,
    const dictionary& dict
)
:
    name_(name),
    patch_(patch),
    traction_
    (
        IOobject
        (
            "traction",
            patch.boundaryMesh().mesh().time().timeName(),
            patch.boundaryMesh().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        patch.boundaryMesh().mesh(),
        dimensionedVector("zero", dimPressure, vector::zero)
    )
{}


Foam::cohesiveZoneModel::cohesiveZoneModel(const cohesiveZoneModel& czm)
:
    name_(czm.name_),
    patch_(czm.patch_),
    traction_(czm.traction_)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


const Foam::surfaceVectorField& Foam::cohesiveZoneModel::meshTraction() const
{
    return traction_;
}


void Foam::cohesiveZoneModel::updateMeshTraction() const
{
    // Reference to the mesh
    const fvMesh& mesh = this->mesh();

    // Face unit normals
    const surfaceVectorField n(mesh.Sf()/mesh.magSf());

    // Lookup the stress field from the solver: this should be up-to-date
    if (mesh.foundObject<surfaceSymmTensorField>("sigmaf"))
    {
        const surfaceSymmTensorField& sigma =
            mesh.lookupObject<surfaceSymmTensorField>("sigmaf");

        traction_ = n & sigma;
    }
    else if (mesh.foundObject<volSymmTensorField>("sigma"))
    {
        const volSymmTensorField& sigma =
            mesh.lookupObject<volSymmTensorField>("sigma");

        traction_ = n & fvc::interpolate(sigma);
    }
    else
    {
        FatalErrorInFunction
            << "Failed in looking up the stress field from the solver"
            << abort(FatalError);
    }
}


// ************************************************************************* //
