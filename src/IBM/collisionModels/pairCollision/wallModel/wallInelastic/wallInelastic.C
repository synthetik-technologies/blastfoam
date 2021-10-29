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

#include "wallInelastic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wallInelastic, 0);
    addToRunTimeSelectionTable
    (
        wallModel,
        wallInelastic,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallInelastic::wallInelastic
(
    const dictionary& dict
)
:
    wallModel(dict, typeName),
    e_(this->coeffDict().template lookup<scalar>("e")),
    staticFrictionCoeff_
    (
        this->coeffDict().template lookup<scalar>("staticFrictionCoeff")
    ),
    dynamicFrictionCoeff_
    (
        this->coeffDict().template lookup<scalar>("dynamicFrictionCoeff")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallInelastic::~wallInelastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallInelastic::evaluateWall
(
    const polyPatch& patch,
    immersedBoundaryObject& object,
    const wallCollisionData& map
) const
{
    label facei(map.hitPointIndex());
    if (facei < 0)
    {
        return;
    }

    const scalar dT = object.pMesh().time().deltaTValue();
    const tensor invI = object.invMomentOfInertia();
    const vector& x = object.centreOfRotation();

    vector hitPoint(map.hitPoint());
    vector normal(map.normal());
    vector v(map.v());

    scalar vN(v & normal);

    if (vN <= 0)
    {
        return;
    }
    vector r = hitPoint - x;
    scalar invMassSum
    (
        1.0/object.mass()
      + (normal & (invI & ((r ^ normal) ^ r)))
    );
    scalar fn = -(1.0 + e_)*vN/invMassSum/dT;
    object.force()[facei] += fn*normal;

    // Frictional force
    vector tangent(v - vN*normal);
    tangent /= max(mag(tangent), small);

    scalar ft(-(v & tangent)/invMassSum/dT);
    if (mag(ft) < small)
    {
        return;
    }

    // Use the effective normal force (gravity, pressure, etc)
    fn = object.forceEff() & normal;
    vector Ft(Zero);
    if (mag(ft) < fn*staticFrictionCoeff_)
    {
        Ft = ft*tangent;
    }
    else
    {
        Ft = -fn*dynamicFrictionCoeff_*tangent;
    }

    // Apply forcing
    object.force()[facei] += Ft;
}


void Foam::wallInelastic::evaluateExternalForce
(
    const polyPatch& patch,
    immersedBoundaryObject& object,
    const wallCollisionData& map
) const
{
    if (!map.size())
    {
        return;
    }

    const label facei = map.hitPointIndex();
    const scalarField& weights(map.weights());
    const vector x = object.centreOfRotation();
    const scalar dT = object.pMesh().time().deltaTValue();

    vector r(map.hitPoint() - x);
    vector normal(map.normal());

    scalar invMassSum
    (
        1.0/object.mass()
      + (normal & (object.invMomentOfInertia() & ((r ^ normal) ^ r)))
    );

    scalar vN
    (
        max
        (
//             max
            (
                map.v() & normal
//                 (object.faceCentres()[facei] - object.faceCentresOld()[facei])/object.pMesh().time().deltaT0Value() & normal
            ),
            0.0
        )
    );

    // Normal effective force
    scalar FeffN(object.forceEff() & normal);
    scalar Fstop(vN/invMassSum/dT);
    vector Fext
    (
        -max(Fstop + FeffN, 0.0)*normal
    );

    Info<<((object.omega() ^ r) & normal)<<" "<<(map.v() & normal)<<endl;
    object.force()[map.hitPointIndex()] += Fext;
    object.forceEff() += Fext;
    object.momentEff() += r ^ Fext;

    //Need to prevent rotation
//     vector Meff(object.momentEff());
//     vector rxn(r ^ normal);
//     scalar FMextN(Zero);
// //     if (rxn[0] > small)
// //     {
// //         FMextN = max(FMextN, Meff[0]/rxn[0]);
// //         Info<<"Mx"<<FMextN<<" "<<Meff[0]<<endl;
// //     }
// //     if (rxn[1] > small)
// //     {
// //         FMextN = max(FMextN, Meff[1]/rxn[1]);
// //         Info<<"My"<<FMextN<<endl;
// //     }
//     if (mag(rxn[2]) > small)
//     {
//         FMextN = Meff[2]/rxn[2];
//     }
//
//     vector FMext(min(-FMextN, 0)*normal);
//
//     object.force()[map.hitPointIndex()] += FMext;
//     object.momentEff() += r ^ FMext;


//     const tensor invI = object.invMomentOfInertia();
//     forAll(map, i)
//     {
//         label facei = map[i][0];
//         label patchFacei = map[i][1];
//         const vector& p = object.faceCentres()[facei];
//         vector normal
//         (
//              patch.faceAreas()[patchFacei]
//             /mag(patch.faceAreas()[patchFacei])
//         );
//         vector r(p - x);
//
//         scalar feff(min(-object.forceEff() & normal, Zero));
//
// //         scalar vN(object.velocity(p) & normal);
// //         if (vN > 0)
// //         {
// //             scalar invMassSum
// //             (
// //                 1.0/object.mass()
// // //               + (normal & (invI & ((r ^ normal) ^ r)))
// //             );
// //             feff = min(feff, -vN/invMassSum/dT);
// //         }
// /*
//         vector rxnxr((r^normal)^r);
//         rxnxr[0] = stabilise(rxnxr[0], 1e-10);
//         rxnxr[1] = stabilise(rxnxr[1], 1e-10);
//         rxnxr[2] = stabilise(rxnxr[2], 1e-10);
//         vector FMeff
//         (
//             cmptDivide(object.momentEff() ^ r, rxnxr)
//         );
//         feff = min(feff, -FMeff & normal);*/
//
//         scalar fn(object.force()[facei] & normal);
//         object.force()[facei] += feff*normal*weights[i];
//     }
}
// ************************************************************************* //
