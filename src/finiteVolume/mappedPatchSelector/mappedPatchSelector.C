/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "mappedPatchSelector.H"
#include "mappedWallFvPatch.H"
#include "mappedMovingWallFvPatch.H"

#include "mappedPolyPatch.H"
#include "mappedWallPolyPatch.H"
#include "mappedVariableThicknessWallPolyPatch.H"

#include "mappedMovingPolyPatch.H"
#include "mappedMovingWallPolyPatch.H"
#include "mappedMovingVariableThicknessWallPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedPatchSelector, 0);
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedPatchSelector::mappedPatchSelector
(
    const fvPatch& pp
)
:
    mappedPatchPtr_(nullptr),
    mappedMovingPatchPtr_(nullptr)
{
    if (isA<mappedWallFvPatch>(pp))
    {
        mappedPatchPtr_ = &refCast<const mappedPatchBase>(pp.patch());
    }
    else
    {
        mappedMovingPatchPtr_ =
            &refCast<const mappedMovingPatchBase>(pp.patch());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedPatchSelector::~mappedPatchSelector()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyMesh& Foam::mappedPatchSelector::sampleMesh() const
{
    if (mappedPatchPtr_)
    {
        return mappedPatchPtr_->sampleMesh();
    }

    return mappedMovingPatchPtr_->sampleMesh();
}


const Foam::polyPatch& Foam::mappedPatchSelector::samplePolyPatch() const
{
    if (mappedPatchPtr_)
    {
        return mappedPatchPtr_->samplePolyPatch();
    }

    return mappedMovingPatchPtr_->samplePolyPatch();
}


Foam::tmp<Foam::pointField> Foam::mappedPatchSelector::samplePoints() const
{
    if (mappedPatchPtr_)
    {
        return mappedPatchPtr_->samplePoints();
    }

    return mappedMovingPatchPtr_->samplePoints();
}


Foam::pointIndexHit Foam::mappedPatchSelector::facePoint
(
    const polyMesh& mesh,
    const label facei,
    const polyMesh::cellDecomposition decompMode
)
{
    if (mappedPatchPtr_)
    {
        return mappedPatchPtr_->facePoint(mesh, facei, decompMode);
    }

    return mappedMovingPatchPtr_->facePoint(mesh, facei, decompMode);
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //
// Function to clear patches if an update has occurred
void Foam::clearMappedPatches(fvMesh& mesh)
{
    forAll(mesh.boundaryMesh(), patchi)
    {
        if (isA<mappedWallPolyPatch>(mesh.boundary()[patchi]))
        {
            polyBoundaryMesh& pbMesh =
                const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
            refCast<mappedPatchBase>(pbMesh[patchi]).clearOut();
        }
        if (isA<mappedMovingWallFvPatch>(mesh.boundary()[patchi]))
        {
            polyBoundaryMesh& pbMesh =
                const_cast<polyBoundaryMesh&>
                (
                    mesh.boundaryMesh()
                );
            refCast<mappedMovingPatchBase>(pbMesh[patchi]).clearOut();
        }
    }
}


bool Foam::setMappedPatchDisplacement(fvMesh& mesh, const word& name)
{
    if (mesh.foundObject<pointVectorField>(name))
    {
        const pointVectorField& D =
            mesh.lookupObject<pointVectorField>(name);

        forAll(mesh.boundaryMesh(), patchi)
        {
            if (isA<mappedMovingWallFvPatch>(mesh.boundary()[patchi]))
            {
                polyBoundaryMesh& pbMesh =
                    const_cast<polyBoundaryMesh&>
                    (
                        mesh.boundaryMesh()
                    );
                refCast<mappedMovingPatchBase>(pbMesh[patchi]).setOffsets(D);
            }
        }
        return true;
    }
    return false;
}


bool Foam::isAMappedType(const fvPatch& patch)
{
    if (isA<mappedWallPolyPatch>(patch))
    {
        return true;
    }
    if (isA<mappedMovingWallFvPatch>(patch))
    {
        return true;
    }
    return false;
}


void Foam::replaceMappedWithMappedMoving(fvMesh& mesh)
{
    List<polyPatch*> newPolyPatches(mesh.boundaryMesh().size());
    polyBoundaryMesh& bm = const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    fvBoundaryMesh& fvBm = const_cast<fvBoundaryMesh&>(mesh.boundary());
    forAll(bm, patchi)
    {
        bool set = false;
        if (isA<mappedPolyPatch>(mesh.boundaryMesh()[patchi]))
        {
            const mappedPolyPatch& patch =
                refCast<const mappedPolyPatch>(mesh.boundaryMesh()[patchi]);

            bm.set(patchi,
                new mappedMovingPolyPatch
                (
                    patch.name(),
                    patch.size(),
                    patch.start(),
                    patch.index(),
                    patch.sampleRegion(),
                    patch.samplePatch(),
                    mesh.boundaryMesh()
                ));
            set = true;
        }
        else if (isA<mappedWallPolyPatch>(mesh.boundaryMesh()[patchi]))
        {
            const mappedWallPolyPatch& patch =
                refCast<const mappedWallPolyPatch>(mesh.boundaryMesh()[patchi]);

            bm.set(patchi,
                new mappedMovingWallPolyPatch
                (
                    patch.name(),
                    patch.size(),
                    patch.start(),
                    patch.index(),
                    patch.sampleRegion(),
                    patch.samplePatch(),
                    mesh.boundaryMesh()
                ));
            set = true;
        }
        else if
        (
            isA<mappedVariableThicknessWallPolyPatch>
            (
                mesh.boundaryMesh()[patchi]
            )
        )
        {
            const mappedVariableThicknessWallPolyPatch& patch =
                refCast<const mappedVariableThicknessWallPolyPatch>
                (
                    mesh.boundaryMesh()[patchi]
                );

            bm.set(patchi,
                new mappedMovingVariableThicknessWallPolyPatch
                (
                    patch.name(),
                    patch.size(),
                    patch.start(),
                    patch.index(),
                    patch.sampleRegion(),
                    patch.samplePatch(),
                    mesh.boundaryMesh()
                ));
            set = true;
        }

        if (set)
        {
            fvBm.set
            (
                patchi,
                new mappedMovingWallFvPatch
                (
                    mesh.boundaryMesh()[patchi],
                    mesh.boundary()
                )
            );
        }
    }

    mesh.clearOut();
}


// ************************************************************************* //
