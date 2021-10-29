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

#include "wallLocalInelastic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wallLocalInelastic, 0);
    addToRunTimeSelectionTable
    (
        wallModel,
        wallLocalInelastic,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallLocalInelastic::wallLocalInelastic
(
    const dictionary& dict
)
:
    wallModel(dict, typeName),
    e_(this->coeffDict().template lookup<scalar>("e"))
{
//     const polyMesh& mesh = cloud.mesh();
//
//     const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
//
//     patchMap_.setSize(bMesh.size(), -1);
//
//     DynamicList<label> wallPatchIndices;
//
//     forAll(bMesh, patchi)
//     {
//         if (isA<wallPolyPatch>(bMesh[patchi]))
//         {
//             wallPatchIndices.append(bMesh[patchi].index());
//         }
//     }
//
//     label nWallPatches = wallPatchIndices.size();
//
//     Estar_.setSize(nWallPatches);
//     Gstar_.setSize(nWallPatches);
//     alpha_.setSize(nWallPatches);
//     b_.setSize(nWallPatches);
//     mu_.setSize(nWallPatches);
//     cohesionEnergyDensity_.setSize(nWallPatches);
//     cohesion_.setSize(nWallPatches);
//
//     scalar maxEstar = -great;
//
//     forAll(wallPatchIndices, wPI)
//     {
//         const dictionary& patchCoeffDict
//         (
//             this->coeffDict().subDict(bMesh[wallPatchIndices[wPI]].name())
//         );
//
//         patchMap_[wallPatchIndices[wPI]] = wPI;
//
//         scalar nu = patchCoeffDict.template lookup<scalar>("poissonsRatio");
//
//         scalar E = patchCoeffDict.template lookup<scalar>("youngsModulus");
//     }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallLocalInelastic::~wallLocalInelastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallLocalInelastic::evaluateWall
(
    const polyPatch& patch,
    immersedBoundaryObject& object,
    const wallCollisionData& map
) const
{
//     vector Fmax = Zero;
//     label index = -1;
//     forAll(map, i)
//     {
//         label facei = map[i][0];
//         label patchFacei = map[i][1];
//         const vector& ofc = object.faceCentres()[facei];
//         const vector& pfc = patch.faceCentres()[patchFacei];
//         scalar patchMagSf(mag(patch.faceAreas()[patchFacei]));
//         vector patchNormal(patch.faceAreas()[patchFacei]/patchMagSf);
//
//         vector v(object.velocity(ofc));
//
//         vector overlap = ofc - pfc;
//         if ((v & patchNormal) > 0)
//         {
//             vector vNorm = patchNormal*(patchNormal & v);
//
//             vector F =
//                 -2.0*e_*object.mass()*vNorm/object.pMesh().time().deltaTValue();
//             if (mag(F) > mag(Fmax))
//             {
//                 Fmax = F;
//                 index = facei;
//             }
//         }
//     }
//     object.force()[index] = Fmax;
}


void Foam::wallLocalInelastic::evaluateExternalForce
(
    const polyPatch& patch,
    immersedBoundaryObject& object,
    const wallCollisionData& map
) const
{}

// ************************************************************************* //
