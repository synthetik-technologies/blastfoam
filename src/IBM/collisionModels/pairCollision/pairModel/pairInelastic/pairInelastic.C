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

#include "pairInelastic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pairInelastic, 0);
    addToRunTimeSelectionTable
    (
        pairModel,
        pairInelastic,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pairInelastic::pairInelastic
(
    const dictionary& dict
)
:
    pairModel(dict, typeName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pairInelastic::~pairInelastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::pairInelastic::evaluatePair
(
    immersedBoundaryObject& objectA,
    immersedBoundaryObject& objectB,
    const pairCollisionData& map
) const
{
    const label facei = map.hitPointIndexA();
    const label facej = map.hitPointIndexB();
    if (facei < 0 && facej < 0)
    {
        return;
    }

    scalar dT = objectA.pMesh().time().deltaTValue();
    scalar massA = objectA.mass();
    scalar massB = objectB.mass();

    tensor invIA = objectA.invMomentOfInertia();
    tensor invIB = objectB.invMomentOfInertia();

    vector xA = objectA.centre();
    vector xB = objectB.centre();

    vector hitPoint(map.hitPoint());
    vector normal(map.normalA());
    vector vAB(map.vAB());


    vector rA(hitPoint - xA);
    vector rB(hitPoint - xB);

    scalar e = 0.5;

    scalar invMassSum
    (
        1.0/massA + 1.0/massB
      + (normal & (invIA & ((rA ^ normal) ^ rA)))
      + (normal & (invIB & ((rB ^ normal) ^ rB)))
    );

    scalar vN(vAB & normal);
    scalar fn(-(1.0 + e)*vN/invMassSum/dT);
    if (facei >= 0)
    {
        objectA.force()[facei] += fn*normal;
    }
    if (facej >= 0)
    {
        objectB.force()[facej] -= fn*normal;
    }

    // Frictional force
    vector tangent(vAB - vN*normal);
    tangent /= max(mag(tangent), small);

    scalar ft(-(vAB & tangent)/invMassSum/dT);
    if (mag(ft) < small)
    {
        return;
    }

    fn = (objectA.forceEff() - objectB.forceEff()) & normal;
    vector Ft(Zero);
    scalar staticFrictionCoeff_ = 0.5;
    scalar dynamicFrictionCoeff_ = 0.5;
    if (mag(ft) < fn*staticFrictionCoeff_)
    {
        Ft = ft*tangent;
    }
    else
    {
        Ft = -fn*dynamicFrictionCoeff_*tangent;
    }

    // Apply forcing
    if (facei >= 0)
    {
        objectA.force()[facei] += Ft;
    }
    if (facej >= 0)
    {
        objectB.force()[facej] -= Ft;
    }
}

// ************************************************************************* //
