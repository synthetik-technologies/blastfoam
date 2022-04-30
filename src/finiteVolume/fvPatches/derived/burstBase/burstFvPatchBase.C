/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "burstFvPatchBase.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "mapPolyMesh.H"
#include "fvPatchMapper.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(burstFvPatchParent, 0);
    defineTypeNameAndDebug(burstFvPatchBase, 0);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::burstFvPatchParent::makeIntact() const
{
    if (!needToRead_)
    {
        return;
    }
    if (mesh_.foundObject<volScalarField>("intact"))
    {
        return;
    }
    needToRead_ = false;

    IOobject intactIO
    (
        "intact",
        mesh_.time().timeName(mesh_.time().startTime().value()),
        mesh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    );

    if (intactIO.typeHeaderOk<volScalarField>(true))
    {
        intactIO.readOpt() = IOobject::MUST_READ;
        intactPtr_.set
        (
            new volScalarField
            (
                intactIO,
                mesh_
            )
        );
    }
    else
    {
        wordList patchTypes(mesh_.boundaryMesh().types());
        wordList boundaryTypes(patchTypes);
        forAll(mesh_.boundary(), patchi)
        {
            if
            (
                isA<burstFvPatchBase>(mesh_.boundary()[patchi])
            || !polyPatch::constraintType
                (
                    mesh_.boundaryMesh()[patchi].type()
                )
            )
            {
                boundaryTypes[patchi] =
                    calculatedFvPatchField<scalar>::typeName;
            }
        }
        intactPtr_.set
        (
            new volScalarField
            (
                intactIO,
                mesh_,
                0.0,
                boundaryTypes,
                patchTypes
            )
        );
        forAll(mesh_.boundaryMesh(), patchi)
        {
            if
            (
                isA<burstPolyPatchBase>(mesh_.boundaryMesh()[patchi])
             || isA<wallPolyPatch>(mesh_.boundaryMesh()[patchi])
            )
            {
                intactPtr_->boundaryFieldRef()[patchi] == 1.0;
            }
        }
    }
    intactPtr_->instance() = mesh_.time().timeName();
}

Foam::burstFvPatchParent::burstFvPatchParent(const fvMesh& mesh)
:
    regIOobject
    (
        IOobject
        (
            typeName,
            mesh.time().timeName(),
            mesh
        )
    ),
    mesh_(mesh),
    needToRead_
    (
        mesh_.time().timeIndex() == mesh_.time().startTimeIndex()
    )
{}


Foam::burstFvPatchParent::~burstFvPatchParent()
{}


Foam::burstFvPatchParent& Foam::burstFvPatchParent::New(const fvMesh& mesh)
{
    if (!mesh.foundObject<burstFvPatchParent>(typeName))
    {
        burstFvPatchParent* parent = new burstFvPatchParent(mesh);
        parent->store(parent);
    }
    return mesh.lookupObjectRef<burstFvPatchParent>(typeName);
}


void Foam::burstFvPatchParent::update
(
    const label patchi,
    const scalarField& intact
)
{
    if (intactPtr_.valid())
    {
        intactPtr_->boundaryFieldRef()[patchi] = intact;
    }
}


void Foam::burstFvPatchBase::updateDeltas()
{
    const fvMesh& mesh(this->patch_.boundaryMesh().mesh());

    //- Clear the interpolation weights
    const_cast<fvMesh&>(mesh).surfaceInterpolation::movePoints();
}


void Foam::burstFvPatchBase::update()
{
    if (patch_.boundaryMesh().mesh().time().timeIndex() == curTimeIndex_)
    {
        return;
    }
    curTimeIndex_ = patch_.boundaryMesh().mesh().time().timeIndex();
    scalarField intact(this->intact());

    if (burstPolyPatch_.update(patch_, intact))
    {
        // Reset partially intact faces
        forAll(intact, i)
        {
            if (intact[i] > 0.01)
            {
                intact[i] = 1;
            }
        }
        parent_.update(patch_.index(), intact);
        updateDeltas();
    }
}
// ************************************************************************* //
