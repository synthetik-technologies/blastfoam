/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "immersedMeshMapper.H"
#include "wedgePolyPatch.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "mappedPatchBase.H"
#include "treeBoundBox.H"
#include "treeDataFace.H"
#include "immersedBoundaryObjectListSolver.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::immersedMeshMapper::interpolationMethod,
        1
    >::names[] =
    {
        "inverseDistance"
    };
}


const
Foam::NamedEnum<Foam::immersedMeshMapper::interpolationMethod, 1>
    Foam::immersedMeshMapper::interpolationMethodNames_;


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::labelList
Foam::immersedMeshMapper::findNearestFaces(const pointField& points) const
{
    (void)immersedMesh_.tetBasePtIs();

    const polyBoundaryMesh& bm = immersedMesh_.boundaryMesh();

     // All the info for nearest. Construct to miss
    List<mappedPatchBase::nearInfo> nearest(points.size());

    const polyPatch& pp = bm[interfaceIndex_];
    pointField faceCentres(immersedMesh_.C().boundaryField()[interfaceIndex_]);
    if (DPtr_ != nullptr)
    {
        faceCentres += (*DPtr_).boundaryField()[interfaceIndex_];
        Info<<faceCentres.size()<<" "<<(*DPtr_).boundaryField()[interfaceIndex_].size()<<endl;
    }


    if (pp.size() > 0)
    {
        labelList bndFaces(pp.size());
        forAll(bndFaces, i)
        {
            bndFaces[i] =  pp.start() + i;
        }

        treeBoundBox overallBb(pp.points());
        overallBb = overallBb.extend(1e-4);

        const indexedOctree<treeDataFace> boundaryTree
        (
            treeDataFace    // all information needed to search faces
            (
                false,                      // do not cache bb
                immersedMesh_,
                bndFaces                    // patch faces only
            ),
            overallBb,                      // overall search domain
            8,                              // maxLevel
            10,                             // leafsize
            3.0                             // duplicity
        );


        forAll(points, i)
        {
            const point sample = points[i];

            scalar span = boundaryTree.bb().mag();

            pointIndexHit info = boundaryTree.findNearest
            (
                sample,
                Foam::sqr(span)
            );

            if (!info.hit())
            {
                info = boundaryTree.findNearest
                (
                    sample,
                    Foam::sqr(great)
                );
            }

            label facei = boundaryTree.shapes().faceLabels()[info.index()];
            const point& fc = faceCentres[facei - pp.start()];

            mappedPatchBase::nearInfo sampleInfo;

            sampleInfo.first() = pointIndexHit
            (
                true,
                fc,
                facei
            );

            sampleInfo.second().first() = magSqr(fc-sample);
            sampleInfo.second().second() = Pstream::myProcNo();

            nearest[i] = sampleInfo;
        }
    }


    // Find nearest.
    Pstream::listCombineGather(nearest, mappedPatchBase::nearestEqOp());
    Pstream::listCombineScatter(nearest);


    // Extract any local faces to sample
    labelList nearestFaces(nearest.size(), -1);

    forAll(nearest, sampleI)
    {
        if (nearest[sampleI].second().second() == Pstream::myProcNo())
        {
            // Store the face to sample
            nearestFaces[sampleI] = nearest[sampleI].first().index();
        }
    }
    return nearestFaces;
}

void Foam::immersedMeshMapper::calcMapping()
{
    //- Make sure patch index is correct
    interfaceIndex_ = immersedMesh_.boundaryMesh()[patchName_].index();
    const polyPatch& patch = immersedMesh_.boundaryMesh()[interfaceIndex_];

    const pointField& meshFaceCentres
    (
        immersedMesh_.C().boundaryField()[interfaceIndex_]
    );
    const pointField& objectFaceCentres(immersedObjectPtr_->faceCentres());
    labelList nearestFaces(findNearestFaces(objectFaceCentres));
    label startIndex = patch.start();

    if (!globalImmersedMeshFaceMapPtr_)
    {
        globalImmersedMeshFaceMapPtr_ = new globalIndex(patch.size());
    }

    immersedObjectToMeshPoints_ = labelListList(patch.size());
    immersedObjectToMeshWeights_ = List<scalarList>(patch.size());
    immersedMeshToObjectPoints_ = labelListList(immersedObjectPtr_->nFaces());
    immersedMeshToObjectWeights_ =
        List<scalarList>(immersedObjectPtr_->nFaces());

    scalarList sumMeshWeights(globalImmersedMeshFaceMapPtr_->size(), 0.0);
    scalarList sumObjectWeights(objectFaceCentres.size(), 0.0);
    forAll(nearestFaces, i)
    {
        // Local face index
        label faceI = nearestFaces[i] - startIndex;
        if (faceI < 0)
        {
            continue;
        }

        // Global face index
        label facei = globalImmersedMeshFaceMapPtr_->toGlobal(faceI);
        scalar invDist =
            1.0
           /max
            (
                mag
                (
                    immersedObjectPtr_->shape().zeroDir
                    (
                        meshFaceCentres[faceI] - objectFaceCentres[i]
                    )
                ),
                small
            );
        immersedObjectToMeshPoints_[faceI].append(i);
        immersedObjectToMeshWeights_[faceI].append(invDist);
        sumMeshWeights[faceI] += invDist;

        immersedMeshToObjectPoints_[i].append(facei);
        immersedMeshToObjectWeights_[i].append(invDist);
        sumObjectWeights[i] += invDist;

        const labelList& faces =
            immersedMesh_.boundaryMesh()[interfaceIndex_].faceFaces()[faceI];
        forAll(faces, j)
        {
            label faceJ = faces[j];
            label facej = globalImmersedMeshFaceMapPtr_->toGlobal(faceJ);
            invDist =
                1.0
               /max(mag(meshFaceCentres[faceJ] - objectFaceCentres[i]), small);
            immersedObjectToMeshPoints_[faceJ].append(i);
            immersedObjectToMeshWeights_[faceJ].append(invDist);
            sumMeshWeights[faceJ] += invDist;

            immersedMeshToObjectPoints_[i].append(facej);
            immersedMeshToObjectWeights_[i].append(invDist);
            sumObjectWeights[i] += invDist;
        }
    }
    forAll(immersedObjectToMeshWeights_, i)
    {
        forAll(immersedObjectToMeshWeights_[i], j)
        {
            immersedObjectToMeshWeights_[i][j] /= sumMeshWeights[i];
        }
    }
    forAll(immersedMeshToObjectWeights_, i)
    {
        reduce(sumObjectWeights[i], sumOp<scalar>());
        forAll(immersedMeshToObjectWeights_[i], j)
        {
            immersedMeshToObjectWeights_[i][j] /= sumObjectWeights[i];
        }
    }
}

// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::immersedMeshMapper::immersedMeshMapper
(
    const fvMesh& immersedMesh,
    const immersedBoundaryObject& immersedObject,
    const word& patchName
)
:
    BalanceMeshObject
    (
        IOobject::groupName("immersedMeshMapper", immersedObject.name()),
        immersedMesh.time()
    ),
    name_(immersedObject.name()),
    pMesh_(immersedObject.pMesh()),
    immersedMesh_(immersedMesh),
    immersedObjectPtr_
    (
        &pMesh_.lookupObject<immersedBoundaryObjectListSolver>
        (
            "immersedBoundaryObjectListSolver"
        ).objects()[name_]
    ),
    mode_(INVERSE_DISTANCE),
    patchName_(patchName),
    immersedObjectToMeshPoints_(),
    immersedObjectToMeshWeights_(),
    immersedMeshToObjectPoints_(immersedObjectPtr_->nFaces()),
    immersedMeshToObjectWeights_(immersedObjectPtr_->nFaces()),
    globalImmersedMeshFaceMapPtr_(nullptr),
    pointDPtr_(nullptr),
    DPtr_(nullptr)
{
    label nFaces = 0;
    wordList patchNames;
    forAll(immersedMesh_.boundaryMesh(), patchi)
    {
        const polyPatch& patch = immersedMesh_.boundaryMesh()[patchi];
        if (patch.name() == patchName)
        {
            interfaceIndex_ = patchi;
            nFaces = patch.size();
            patchNames.clear();
            patchNames.resize(1, patch.name());
            break;
        }
        else if
        (
            !isA<wedgePolyPatch>(patch)
         && !isA<emptyPolyPatch>(patch)
         && !isA<processorPolyPatch>(patch)
        )
        {
            interfaceIndex_ = patchi;
            nFaces += patch.size();
            patchNames.append(patch.name());
            if (patchName_ == word::null)
            {
                patchName_ = patch.name();
            }
        }
    }
    if (patchNames.size() > 1)
    {
        FatalErrorInFunction
            << "Immersed meshes can only have 1 interface patch" << nl
            << patchNames.size() << " have been specified" << nl
            << patchNames << nl
            << abort(FatalError);
    }
    immersedObjectToMeshPoints_.resize(nFaces);
    immersedObjectToMeshWeights_.resize(nFaces);

    calcMapping();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::immersedMeshMapper::~immersedMeshMapper()
{
    deleteDemandDrivenData(globalImmersedMeshFaceMapPtr_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::immersedMeshMapper::setDisplacement
(
    const pointVectorField& pointD,
    const volVectorField& D
)
{
    pointDPtr_ = &pointD;
    DPtr_ = &D;
}

// ************************************************************************* //
