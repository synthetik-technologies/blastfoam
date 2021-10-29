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

#include "pairCollision.H"
#include "pairModel.H"
#include "wallModel.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pairCollision, 0);
    addToRunTimeSelectionTable
    (
        collisionModel,
        pairCollision,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pairCollision::objectInteraction()
{
    forAll(objects_, i)
    {
        for (label j = i+1; j < objects_.size(); j++)
        {
            Info<<objects_[i].name()<<" "<<objects_[j].name()<<endl;
            pairCollisionData map(objects_[i], objects_[j]);
            map.calcMapping();

            pairModel_->evaluatePair
            (
                objects_[i],
                objects_[j],
                map
            );
        }
    }
}


void Foam::pairCollision::wallInteraction()
{
    const polyMesh& mesh = objects_[0].pMesh();

    forAll(patchNames_, patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchNames_[patchi]];
        boundBox patchBb
        (
            IndirectList<point>(mesh.points(), patch.meshPoints())()
        );
        vector avgNormal
        (
            gSum(patch.faceAreas())
           /gSum(mag(patch.faceAreas()))
        );
        patchBb.min() += min(avgNormal, vector::zero)*great;
        patchBb.max() += max(avgNormal, vector::zero)*great;

        forAll(objects_, objecti)
        {
            wallObjectMaps_[patchi][objecti].patchBB() = patchBb;
            wallObjectMaps_[patchi][objecti].calcMapping();

            wallModel_->evaluateWall
            (
                patch,
                objects_[objecti],
                wallObjectMaps_[patchi][objecti]
            );
        }

//         // Update external forces and moments
//         forAll(objects_, objecti)
//         {
//             //- Save old forces and moments and set to zero to get the
//             //  forces and moments added from the wall
//             vector Feff = returnReduce
//             (
//                 objects_[objecti].forceEff(),
//                 maxMagSqrOp<vector>()
//             );
//
//             // Use the moment from the processor with the largest
//             // force
//             if (mag(Feff - objects_[objecti].forceEff()) > small)
//             {
//                 objects_[objecti].momentEff() = Zero;
//             }
//             reduce(objects_[objecti].momentEff(), sumOp<vector>());
//             objects_[objecti].forceEff() = Feff + forceEffs[objecti];
//             objects_[objecti].momentEff() += momentEffs[objecti];
//
//
//             //- Save old forces and moments and set to zero to get the
//             //  forces and moments added from the wall
//             vector Fext = returnReduce
//             (
//                 objects_[objecti].forceExt(),
//                 maxMagSqrOp<vector>()
//             );
//
//             // Use the moment from the processor with the largest
//             // force
//             if (mag(Fext - objects_[objecti].forceExt()) > small)
//             {
//                 objects_[objecti].momentExt() = Zero;
//             }
//             reduce(objects_[objecti].momentExt(), sumOp<vector>());
//             objects_[objecti].forceExt() = Fext + forceExts[objecti];
//             objects_[objecti].momentExt() += momentExts[objecti];
//         }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pairCollision::pairCollision
(
    const polyMesh& mesh,
    PtrList<immersedBoundaryObject>& objects,
    const dictionary& dict
)
:
    collisionModel(mesh, objects, dict),
    pairModel_(pairModel::New(dict)),
    wallModel_(wallModel::New(dict)),
    wallObjectMaps_()
{

    forAll(mesh.boundaryMesh(), patchi)
    {
        if (isA<wallPolyPatch>(mesh.boundaryMesh()[patchi]))
        {
            patchNames_.append(mesh.boundaryMesh()[patchi].name());
        }
    }
    wallObjectMaps_.resize(patchNames_.size());
    forAll(patchNames_, i)
    {
        wallObjectMaps_.set
        (
            i,
            new PtrList<wallCollisionData>(objects.size())
        );
        forAll(objects, j)
        {
            wallObjectMaps_[i].set
            (
                j,
                new wallCollisionData
                (
                    mesh.boundaryMesh()[patchNames_[i]],
                    objects_[j]
                )
            );
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pairCollision::~pairCollision()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pairCollision::collide()
{
    objectInteraction();
    wallInteraction();
}


void Foam::pairCollision::updateExternalForce()
{
    const polyMesh& mesh = objects_[0].pMesh();

    forAll(patchNames_, patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchNames_[patchi]];
        forAll(objects_, objecti)
        {
            wallModel_->evaluateExternalForce
            (
                patch,
                objects_[objecti],
                wallObjectMaps_[patchi][objecti]
            );
        }
    }
}


// ************************************************************************* //
