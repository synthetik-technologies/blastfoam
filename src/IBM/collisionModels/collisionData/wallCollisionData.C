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

#include "wallCollisionData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallCollisionData::wallCollisionData
(
    const polyPatch& patch,
    const immersedBoundaryObject& object
)
:
    labelPairList(0),
    patch_(patch),
    object_(object),
    weights_(0),
    hitPoint_(Zero),
    normal_(Zero),
    v_(Zero),
    hitPointIndex_(-1),
    reducedMap_(0)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallCollisionData::~wallCollisionData()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallCollisionData::calcMapping()
{
    this->clear();
    labelPairList& map(*this);

    // Check if bounding boxes intersect
    const boundBox& objectBb(object_.bounds());
    if (objectBb.overlaps(patchBB_))
    {
        const pointField& objectPoints(object_.faceCentres());

        // Create list of points in objectBb
        labelList validPatchPoints
        (
            object_.calcInside(patch_.faceCentres())
        );
        if (validPatchPoints.size() > 0)
        {
            label pj = 0;
            labelList validObjectPoints(objectPoints.size());
            forAll(objectPoints, pointj)
            {
                if (patchBB_.contains(objectPoints[pointj]))
                {
                    validObjectPoints[pj++] = pointj;
                }
            }

            if (pj > 0)
            {
                validObjectPoints.resize(pj);

                label minIndex = 0;
                scalar minDist = great;
                scalar dist = 0;

                // Find nearest faces for points in object j
                map.resize(validObjectPoints.size());
                forAll(validObjectPoints, pointj)
                {
                    minIndex = 0;
                    minDist = great;
                    forAll(validPatchPoints, pointi)
                    {
                        dist =
                            mag
                            (
                                patch_.faceCentres()[validPatchPoints[pointi]]
                            - objectPoints[validObjectPoints[pointj]]
                            );
                        if (dist < minDist)
                        {
                            minDist = dist;
                            minIndex = validPatchPoints[pointi];
                        }
                    }
                    map[pointj][0] = validObjectPoints[pointj];
                    map[pointj][1] = minIndex;
                }
            }
        }
    }

    reducedWeights_.resize(map.size());
    reducedWeights_ = 0.0;
    weights_.resize(map.size());
    weights_ = 0.0;
    hitPoint_ = Zero;
    normal_ = Zero;
    v_ = Zero;
    hitPointIndex_ = -1;

    if (returnReduce(map.size(), sumOp<label>()) == 0)
    {
        return;
    }

    scalar area = 0.0;

    scalar minR = great;
    label wi = 0;
    reducedMap_.resize(map.size());
    forAll(map, i)
    {
        label facei = map[i][0];
        label patchFacei = map[i][1];
        const vector& p = object_.faceCentres()[facei];
        const vector& fc = patch_.faceCentres()[patchFacei];
        const scalar magSf(mag(object_.Sf()[facei]));

        vector normal
        (
            patch_.faceAreas()[patchFacei]
           /mag(patch_.faceAreas()[patchFacei])
        );
        vector v(object_.v(p));

         area += magSf;
        hitPoint_ += magSf*p;
        normal_ += normal*magSf;

        scalar r(mag(fc - object_.centre()));
        if (r < minR)
        {
            minR = r;
            hitPointIndex_ = facei;
        }

        scalar overlap(mag(p - fc));
        weights_[i] = overlap;

        if ((v & normal) > 0)
        {
            reducedWeights_[wi] = overlap;
            reducedMap_[wi++] = map[i];
        }
    }

    reduce(area, sumOp<scalar>());
    if (area < small)
    {
        return;
    }

    scalar test(minR);
    label proc(Pstream::myProcNo());
    reduce(test, minOp<scalar>());
    if (minR != test)
    {
        proc = -1;
    }
    reduce(proc, maxOp<label>());
    if (proc != Pstream::myProcNo())
    {
        hitPointIndex_ = -1;
    }

    reduce(hitPoint_, sumOp<vector>());
    reduce(normal_, sumOp<vector>());
    reduce(v_, sumOp<vector>());

    hitPoint_ /= area;
    normal_ /= area;
    v_ = object_.v(hitPoint_);

    // Return the true weights across all processors
    reducedWeights_.resize(wi);
    reducedMap_.resize(wi);
    reducedWeights_ /=
        returnReduce(sum(reducedWeights_), sumOp<scalar>());

    weights_ /= returnReduce(sum(weights_), sumOp<scalar>());
}

// ************************************************************************* //
