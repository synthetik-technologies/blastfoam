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

#include "pairCollisionData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pairCollisionData::pairCollisionData
(
    const immersedBoundaryObject& objectA,
    const immersedBoundaryObject& objectB
)
:
    labelPairList(0),
    objectA_(objectA),
    objectB_(objectB),
    weights_(0),
    reducedWeights_(0),
    hitPoint_(Zero),
    normalA_(Zero),
    normalB_(Zero),
    vAB_(Zero),
    hitPointIndexA_(-1),
    hitPointIndexB_(-1),
    reducedMap_(0)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pairCollisionData::~pairCollisionData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pairCollisionData::calcMapping()
{
    this->clear();
    labelPairList& map(*this); // closest faces A&B

    // Check if bounding boxes intersect
    if (objectA_.bounds().overlaps(objectB_.bounds()))
    {
        //- Only use points inside the current mesh
        const pointField& pointsA(objectA_.faceCentres());

        // Create list of points in objectB
        labelList validPointsA(objectB_.calcInside(pointsA));
        if (validPointsA.size() > 0)
        {
            const pointField& pointsB(objectB_.faceCentres());

            // Create list of points in objectA
            labelList validPointsB(objectA_.calcInside(pointsB));
            if (validPointsB.size() > 0)
            {
                // Find nearest faces for points in object A
                map.resize(validPointsA.size());

                label minIndex = 0;
                scalar minDist = great;
                scalar dist = 0;

                forAll(validPointsA, pointi)
                {
                    minIndex = -1;
                    minDist = great;
                    forAll(validPointsB, pointj)
                    {
                        dist =
                            mag
                            (
                                pointsA[validPointsA[pointi]]
                              - pointsB[validPointsB[pointj]]
                            );
                        if (dist < minDist)
                        {
                            minDist = dist;
                            minIndex = validPointsB[pointj];
                        }
                    }
                    map[pointi][0] = validPointsA[pointi];
                    map[pointi][1] = minIndex;
                }
            }
        }
    }

    reducedWeights_.resize(map.size());
    weights_ = 0.0;
    reducedWeights_ = 0.0;
    hitPoint_ = Zero;
    normalA_ = Zero;
    normalB_ = Zero;
    vAB_ = Zero;
    hitPointIndexA_ = -1;
    hitPointIndexB_ = -1;

    if (returnReduce(map.size(), sumOp<label>()) == 0)
    {
        return;
    }

    scalar areaAB = 0.0;
    scalar areaA = 0.0;
    scalar areaB = 0.0;

    scalar minRA = great;
    scalar minRB = great;
    label wi = 0;
    reducedMap_.resize(map.size());
    forAll(map, i)
    {
        label facei = map[i][0];
        label facej = map[i][1];
        const vector& pA = objectA_.faceCentres()[facei];
        const vector& pB = objectB_.faceCentres()[facej];

        vector SfA(objectA_.Sf()[facei]);
        scalar magSfA(mag(SfA));
        vector SfB(objectB_.Sf()[facej]);
        scalar magSfB(mag(SfB));
        scalar magSfAB(magSfA*magSfB);

        areaA += magSfA;
        areaB += magSfB;
        areaAB += magSfAB;

        hitPoint_ += magSfAB*0.5*(pA + pB);

        normalA_ += SfA;
        normalB_ += SfB;

        scalar rA(mag(pB - objectA_.centre()));
        if (rA < minRA)
        {
            minRA = rA;
            hitPointIndexA_ = facei;
        }

        scalar rB(mag(pA - objectB_.centre()));
        if (rB < minRB)
        {
            minRB = rB;
            hitPointIndexB_ = facej;
        }

        vector vAB(objectA_.v(pA) - objectB_.v(pB));

        if ((vAB & SfA) < 0 || (vAB & SfB) > 0)
        {
            reducedWeights_[wi] = mag(pA - pB);
            reducedMap_[wi++] = map[i];
        }
    }

    reduce(areaAB, sumOp<scalar>());
    if (areaAB < small)
    {
        return;
    }

    //- Get face index for objectA
    scalar test(minRA);
    label proc(Pstream::myProcNo());
    reduce(test, minOp<scalar>());
    if (minRA != test)
    {
        proc = -1;
    }
    reduce(proc, maxOp<label>());
    if (proc != Pstream::myProcNo())
    {
        hitPointIndexA_ = -1;
    }

    //- Get face index for objectB
    test = minRB;
    proc = Pstream::myProcNo();
    reduce(test, minOp<scalar>());
    if (minRB != test)
    {
        proc = -1;
    }
    reduce(proc, maxOp<label>());
    if (proc != Pstream::myProcNo())
    {
        hitPointIndexB_ = -1;
    }


    reduce(hitPoint_, sumOp<vector>());
    reduce(normalA_, sumOp<vector>());
    reduce(normalB_, sumOp<vector>());
    reduce(vAB_, sumOp<vector>());
    reduce(areaA, sumOp<scalar>());
    reduce(areaB, sumOp<scalar>());

    hitPoint_ /= areaAB;
    normalA_ /= areaA;
    normalB_ /= areaB;
    vAB_ = objectA_.v(hitPoint_) - objectB_.v(hitPoint_);

    // Return the true weights across all processors
    reducedWeights_.resize(wi);
    reducedMap_.resize(wi);
    weights_ /= returnReduce(sum(weights_), sumOp<scalar>());
    reducedWeights_ /=
        returnReduce(sum(reducedWeights_), sumOp<scalar>());
}

// ************************************************************************* //
