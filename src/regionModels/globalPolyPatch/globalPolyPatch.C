/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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

\*---------------------------------------------------------------------------*/

#include "globalPolyPatch.H"
#include "coupledGlobalPolyPatch.H"
#include "polyPatchID.H"
#include "volFields.H"
#include "pointPatchFields.H"
#include "valuePointPatchFields.H"
#include "globalPoints.H"
#include "vtkWritePolyData.H"
#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(globalPolyPatch, 0);
    defineRunTimeSelectionTable(globalPolyPatch, patch);
    addToRunTimeSelectionTable
    (
        globalPolyPatch,
        globalPolyPatch,
        patch
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::globalPolyPatch::calcPhysicalPatch() const
{
    if (debug)
    {
        InfoInFunction
            << "Calculating physical patch"
            << endl;
    }

    if (physicalPatchPtr_.valid())
    {
        FatalErrorInFunction
            << "Physical patch already calculated"
            << abort(FatalError);
    }

    // Patch index
    const label patchID = polyPatch_.index();

    // Allocate old patch points
    if (mesh_.moving() || displacementField_ != "none")
    {
        displacedPoints0Ptr_.reset
        (
            new pointField(mesh_.oldPoints(), polyPatch_.meshPoints())
        );
    }

    // Insert my points
    if (displacementField_ != "none")
    {
        displacedPointsPtr_.reset(new pointField(polyPatch_.localPoints()));
        tmp<vectorField> tppointD;
        tmp<vectorField> tppointD0;
        if (this->mesh_.foundObject<volVectorField>(displacementField_))
        {
            const volVectorField& D =
                this->mesh_.lookupObject<volVectorField>
                (
                    displacementField_
                );

            tppointD = faceToPoint(D.boundaryField()[patchID]);
            if (displacedPoints0Ptr_.valid())
            {
                tppointD0 =
                    faceToPoint
                    (
                        D.oldTime().boundaryField()[patchID]
                    );
            }
        }
        else if
        (
            this->mesh_.foundObject<pointVectorField>(displacementField_)
        )
        {
            const pointVectorField& pointD =
                this->mesh_.lookupObject<pointVectorField>
                (
                    displacementField_
                );
            const pointPatchVectorField& ppointD =
                pointD.boundaryField()[patchID];
            if (isA<valuePointPatchVectorField>(ppointD))
            {
                tppointD = tmp<vectorField>
                (
                    dynamicCast<const valuePointPatchVectorField>(ppointD)
                );
                if (displacedPoints0Ptr_.valid())
                {
                    tppointD0 = tmp<vectorField>
                    (
                        dynamicCast<const valuePointPatchVectorField>
                        (
                            pointD.oldTime().boundaryField()[patchID]
                        )
                    );
                }
            }
            else
            {
                tppointD = ppointD.patchInternalField();
                if (displacedPoints0Ptr_.valid())
                {
                    tppointD0 =
                        pointD.oldTime().boundaryField()
                        [
                            patchID
                        ].patchInternalField();
                }
            }
        }
        else
        {
            FatalErrorInFunction
                << "Could not find " << displacementField_
                << "in " << mesh_.name() << " region." << endl
                << abort(FatalError);
        }

        if (inverseDisplacement_)
        {
            displacedPointsPtr_() -= tppointD;
            if (tppointD0.valid())
            {
                displacedPoints0Ptr_() -= tppointD0;
            }

        }
        else
        {
            displacedPointsPtr_() += tppointD;
            if (tppointD0.valid())
            {
                displacedPoints0Ptr_() += tppointD0;
            }
        }
    }

    physicalPatchPtr_.reset
    (
        new primitivePatch
        (
            SubList<face>(polyPatch_.localFaces(), polyPatch_.size()),
            displacedPointsPtr_.valid()
          ? displacedPointsPtr_()
          : polyPatch_.localPoints()
        )
    );

    if (mesh_.moving())
    {
        physicalPatch0Ptr_.reset
        (
            new primitivePatch
            (
                SubList<face>(polyPatch_.localFaces(), polyPatch_.size()),
                displacedPoints0Ptr_()
            )
        );
    }
}


void Foam::globalPolyPatch::calcPointToFaceInterpolation() const
{
    if (debug)
    {
        InfoInFunction
            << "Calculating point to face interpolation weights"
            << endl;
    }

    if (pointToFaceInterpolatorPtr_.valid())
    {
        FatalErrorInFunction
            << "Face to point weights already set"
            << abort(FatalError);
    }

    pointToFaceInterpolatorPtr_.reset
    (
        new primitivePatchInterpolation(physicalPatch())
    );
}


void Foam::globalPolyPatch::calcFaceToPointInterpolation() const
{
    if (debug)
    {
        InfoInFunction
            << "Calculating face to patch interpolation weights"
            << endl;
    }

    if (faceToPointWeightsPtr_ || faceToPointSumWeightsPtr_)
    {
        FatalErrorInFunction
            << "Face to point weights already set"
            << abort(FatalError);
    }

    const primitivePatch& patch = physicalPatch();
    const labelListList& pointFaces = patch.pointFaces();
    const pointField& points = patch.points();
    const pointField& faceCentres = patch.faceCentres();

    // Compute face to point weights (inverse distance)
    faceToPointWeightsPtr_ = new List<scalarField>(points.size());
    List<scalarField>& weights = *faceToPointWeightsPtr_;

    faceToPointSumWeightsPtr_ = new scalarField(points.size(), 0.0);
    scalarField& sumWeights = *faceToPointSumWeightsPtr_;

    forAll(points, pointi)
    {
        const labelList& pfs = pointFaces[pointi];
        scalarField& ws = weights[pointi];
        ws.setSize(pfs.size());
        scalar sumW = 0.0;
        forAll(pfs, pfi)
        {
            const label facei = pfs[pfi];
            const scalar w = 1.0/mag(points[pointi] - faceCentres[facei]);
            ws[pfi] = w;
            sumW += w;
        }
        sumWeights[pointi] = sumW;
    }

    syncTools::syncPointList
    (
        mesh_,
        polyPatch_.meshPoints(),
        sumWeights,
        plusEqOp<scalar>(),
        0.0
    );
}


const Foam::List<Foam::scalarField>&
Foam::globalPolyPatch::faceToPointWeights() const
{
    if (!faceToPointWeightsPtr_)
    {
        calcFaceToPointInterpolation();
    }
    return *faceToPointWeightsPtr_;
}



const Foam::scalarField&
Foam::globalPolyPatch::faceToPointSumWeights() const
{
    if (!faceToPointSumWeightsPtr_)
    {
        calcFaceToPointInterpolation();
    }
    return *faceToPointSumWeightsPtr_;
}


void Foam::globalPolyPatch::clearOut() const
{
    displacedPointsPtr_.clear();
    physicalPatchPtr_.clear();
    displacedPoints0Ptr_.clear();
    physicalPatch0Ptr_.clear();

    pointToFaceInterpolatorPtr_.clear();

    deleteDemandDrivenData(faceToPointWeightsPtr_);
    deleteDemandDrivenData(faceToPointSumWeightsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::globalPolyPatch::globalPolyPatch
(
    const dictionary& dict,
    const polyPatch& patch
)
:
    mesh_(patch.boundaryMesh().mesh()),
    patchName_(patch.name()),
    polyPatch_(patch),
    displacementField_
    (
        dict.lookupOrDefault<word>("displacementField", "none")
    ),
    inverseDisplacement_(false),

    displacedPointsPtr_(),
    physicalPatchPtr_(),
    displacedPoints0Ptr_(),
    physicalPatch0Ptr_(),
    pointToFaceInterpolatorPtr_(),
    faceToPointWeightsPtr_(nullptr),
    faceToPointSumWeightsPtr_(nullptr)
{}


Foam::globalPolyPatch::globalPolyPatch
(
    const polyPatch& patch,
    const word& displacementField
)
:
    mesh_(patch.boundaryMesh().mesh()),
    patchName_(patch.name()),
    polyPatch_(mesh_.boundaryMesh()[mesh_.boundaryMesh().findPatchID(patchName_)]),
    displacementField_(displacementField),
    inverseDisplacement_(false),
    displacedPointsPtr_(),
    physicalPatchPtr_(),
    displacedPoints0Ptr_(),
    physicalPatch0Ptr_(),
    pointToFaceInterpolatorPtr_(),
    faceToPointWeightsPtr_(nullptr),
    faceToPointSumWeightsPtr_(nullptr)
{}


Foam::autoPtr<Foam::globalPolyPatch> Foam::globalPolyPatch::New
(
    const dictionary& dict,
    const polyPatch& patch
)
{
    patchConstructorTable::iterator cstrIter =
        patchConstructorTablePtr_->find(patch.type());

    if (cstrIter != patchConstructorTablePtr_->end())
    {
        return cstrIter()(dict, patch);
    }
    return autoPtr<globalPolyPatch>
    (
        new globalPolyPatch(dict, patch)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::globalPolyPatch::~globalPolyPatch()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyMesh& Foam::globalPolyPatch::mesh() const
{
    return mesh_;
}


const Foam::primitivePatch& Foam::globalPolyPatch::physicalPatch() const
{
    if (!physicalPatchPtr_.valid())
    {
        calcPhysicalPatch();
    }

    return physicalPatchPtr_();
}


const Foam::primitivePatch& Foam::globalPolyPatch::physicalPatch0() const
{
    if (!physicalPatchPtr_.valid())
    {
        calcPhysicalPatch();
    }

    return
        physicalPatch0Ptr_.valid()
      ? physicalPatch0Ptr_()
      : physicalPatchPtr_();
}


const Foam::primitivePatchInterpolation&
Foam::globalPolyPatch::pointToFaceInterpolator() const
{
    if (!pointToFaceInterpolatorPtr_.valid())
    {
        if (debug)
        {
            InfoInFunction
                << "Calculating local patch interpolator"
                << endl;
        }
        pointToFaceInterpolatorPtr_.reset
        (
            new primitivePatchInterpolation(physicalPatch())
        );
    }

    return pointToFaceInterpolatorPtr_();
}


void Foam::globalPolyPatch::update()
{
    physicalPatch();
}


void Foam::globalPolyPatch::updateMesh()
{
    clearOut();
}


void Foam::globalPolyPatch::updateMeshOther() const
{}


void Foam::globalPolyPatch::movePoints(const bool clear)
{
    if (clear)
    {
        clearOut();
    }
}


void Foam::globalPolyPatch::movePointsOther(const bool clear) const
{}


bool Foam::globalPolyPatch::write() const
{
    bool good = true;
//     if (debug && mesh_.time().outputTime())
//     {
//         mkDir("VTK");
//         globalPatch().writeVTK
//         (
//             "VTK/"
//             + patch_.name() + '_'
//             + Foam::name(mesh_.time().timeIndex())
//         );
//     }
    return returnReduce(good, andOp<bool>());
}

// ************************************************************************* //
