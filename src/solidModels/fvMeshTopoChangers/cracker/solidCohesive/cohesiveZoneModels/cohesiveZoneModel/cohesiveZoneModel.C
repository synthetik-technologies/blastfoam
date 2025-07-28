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
    defineTypeNameAndDebug(cohesiveZoneModelMaster, 0);
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


Foam::cohesiveZoneModelMaster::cohesiveZoneModelMaster
(
    const fvPatch& p
)
:
    patch_(p),
    dict_(NULL)
{}


Foam::cohesiveZoneModelMaster::cohesiveZoneModelMaster
(
    const fvPatch& p,
    const dictionary& dict
)
:
    patch_(p),
    dict_(dict)
{}


Foam::cohesiveZoneModelMaster::cohesiveZoneModelMaster
(
    const fvPatch& p,
    const cohesiveZoneModelMaster& czm
)
:
    patch_(p),
    dict_(czm.dict_)
{
    if (czm.cohesiveZoneModelPtr_.valid())
    {
        cohesiveZoneModelPtr_.set
        (
            czm.cohesiveZoneModelPtr_->clone().ptr()
        );
    }
}


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


void Foam::cohesiveZoneModelMaster::calcCohesiveZone() const
{
    if (cohesiveZoneModelPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    cohesiveZoneModelPtr_ =
        cohesiveZoneModel::New
        (
            "type", patch_, dict_.subDict("cohesiveZoneModel")
        );
}


const Foam::cohesiveZoneModel&
Foam::cohesiveZoneModelMaster::cohesiveZone() const
{
    if (!cohesiveZoneModelPtr_.valid())
    {
        calcCohesiveZone();
    }

    return cohesiveZoneModelPtr_();
}


Foam::cohesiveZoneModel&
Foam::cohesiveZoneModelMaster::cohesiveZone()
{
    if (!cohesiveZoneModelPtr_.valid())
    {
        calcCohesiveZone();
    }

    return cohesiveZoneModelPtr_();
}

void Foam::cohesiveZoneModelMaster::autoMap
(
    const fvPatchFieldMapper& m
)
{
    if (cohesiveZoneModelPtr_.valid())
    {
        cohesiveZoneModelPtr_->autoMap(m);
    }
}


void Foam::cohesiveZoneModelMaster::rmap
(
    const cohesiveZoneModelMaster& czmm,
    const labelList& addr
)
{
    if (cohesiveZoneModelPtr_.valid())
    {
        cohesiveZoneModelPtr_->rmap(czmm.cohesiveZone(), addr);
    }
}

void Foam::cohesiveZoneModelMaster::write(Ostream& os) const
{
    if (cohesiveZoneModelPtr_.valid())
    {
        os.writeKeyword("cohesiveZoneModel") << nl;
        os << indent << token::BEGIN_BLOCK << nl;
        cohesiveZone().write(os);
        os << indent << token::END_BLOCK << nl;
    }
}

// ************************************************************************* //
