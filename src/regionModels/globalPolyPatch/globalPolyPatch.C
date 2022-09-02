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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(globalPolyPatch, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::globalPolyPatch::calcGlobalPatch() const
{
    if (debug)
    {
        InfoInFunction
            << "Calculating primitive patch"
            << endl;
    }

    if
    (
        globalPatchPtr_.valid()
     || pointToGlobalAddrPtr_.valid()
     || faceToGlobalAddrPtr_.valid()
    )
    {
        FatalErrorInFunction
            << "primitive face zone patch and addressing already calculated"
            << abort(FatalError);
    }

    // Get patch
    polyPatchID patchID
    (
        patchName_,
        mesh_.boundaryMesh()
    );
    const polyPatch& patch = mesh_.boundaryMesh()[patchID.index()];

    if (!patchID.active())
    {
        FatalErrorInFunction
            << "Cannot find patch " << patchName_
            << abort(FatalError);
    }

    // Collect points and faces from all processors
    typedef List<point> pointList;
    typedef List<pointList> pointListList;

    pointListList procPatchPoints(Pstream::nProcs());
    faceListList procPatchFaces(Pstream::nProcs());
    labelListList ownerPoint(Pstream::nProcs());

    globalPoints gPoints(mesh_, true, false);
    const globalIndex gIndex(mesh_.nPoints());
    const Map<label>& meshToProcPoint = gPoints.meshToProcPoint();
    const DynamicList<labelPairList>& procPoints = gPoints.procPoints();

    // Add points and faces if the patch is not empty
    if (!mesh_.boundaryMesh()[patchID.index()].empty())
    {
        ownerPoint[Pstream::myProcNo()].setSize(patch.nPoints(), -1);
        const labelList& meshPoints = patch.meshPoints();
        forAll(meshPoints, pi)
        {
            const label pointi = meshPoints[pi];
            Map<label>::const_iterator iter = meshToProcPoint.find(pointi);
            if (iter != meshToProcPoint.cend())
            {
                ownerPoint[Pstream::myProcNo()][pi] =
                    gIndex.toGlobal
                    (
                        procPoints[iter()][0].second(),
                        procPoints[iter()][0].first()
                    );
            }
            else
            {
                ownerPoint[Pstream::myProcNo()][pi] =
                    gIndex.toGlobal(pointi);
            }
        }

        // Insert my points
        pointField pts(patch.localPoints());
        if (displacementField_ != "none")
        {
            if (this->mesh_.foundObject<volVectorField>(displacementField_))
            {
                PrimitivePatchInterpolation<polyPatch> patchInterp
                (
                    this->patch_
                );
                pts +=
                    patchInterp.faceToPointInterpolate
                    (
                        this->mesh_.lookupObject<volVectorField>
                        (
                            displacementField_
                        ).boundaryField()[this->patch_.index()]
                    );
            }
            else if
            (
                this->mesh_.foundObject<pointVectorField>(displacementField_)
            )
            {
                const pointPatchVectorField& disp =
                    this->mesh_.lookupObject<pointVectorField>
                    (
                        displacementField_
                    ).boundaryField()[this->patch_.index()];
                if (isA<valuePointPatchVectorField>(disp))
                {
                    pts +=
                        dynamicCast<const valuePointPatchVectorField>(disp);
                }
                else
                {
                    pts += disp.patchInternalField();
                }
            }
            else
            {
                FatalErrorInFunction
                    << "Could not find " << displacementField_
                    << "in " << mesh_.name() << " region." << endl
                    << abort(FatalError);
            }
        }

        // Insert my points
        procPatchPoints[Pstream::myProcNo()].transfer(pts);

        // Insert my faces
        procPatchFaces[Pstream::myProcNo()] =
            mesh_.boundaryMesh()[patchID.index()].localFaces();
    }

    // Communicate points
    Pstream::gatherList(procPatchPoints);
    Pstream::scatterList(procPatchPoints);

    // Communicate faces
    Pstream::gatherList(procPatchFaces);
    Pstream::scatterList(procPatchFaces);

    // Communicate point owners
    Pstream::gatherList(ownerPoint);
    Pstream::scatterList(ownerPoint);

    // At this point, all points and faces for the current patch
    // are available.

    // Count the number of faces in the global face zone
    label nZoneFaces = 0;
    label nZonePoints = 0;

    // Sum up points and faces to add
    forAll (procPatchFaces, procI)
    {
        nZonePoints += procPatchPoints[procI].size();
        nZoneFaces += procPatchFaces[procI].size();
    }

    if (debug)
    {
        Info<< "Global zone size for patch " << patchID.name()
            << ": " << nZoneFaces << endl;
    }

    if (nZoneFaces == 0)
    {
        FatalErrorInFunction
            << "Patch " << patchID.name()
            << " appears to be globally empty.  "
            << "Please check definition."
            << abort(FatalError);
    }

    // Record current points and faces to add
    pointField zonePoints(nZonePoints);
    faceList zoneFaces(nZoneFaces);

    // PC, 15/12/17
    // I will keep track of duplicate points with this set
    label nDuplicatePoints = 0;

    label nCurPoints = 0;
    label nCurFaces = 0;

    Map<label> procPointMap(mesh_.nPoints());

    // Collect all points and faces
    forAll(procPatchPoints, procI)
    {
        // Add points from all processors except self
        const pointList& curProcPoints = procPatchPoints[procI];

        // Label point map for the current processor
        labelList pointMap(curProcPoints.size());

        // Add points from all processors
        forAll(curProcPoints, pointI)
        {
            // Note: possible removal of duplicate points here
            // HJ, 28/Dec/2016
            // PC, 15/Dec/2017: the edge loops are wrong unless the duplicates
            // are removed, because the internal patch faces will not be
            // connected across the processor boundaries. So I will keep track
            // of the duplicates with a HashTable and remove them

            // Current point
            const point& curPoint = curProcPoints[pointI];

            if (procPointMap.insert(ownerPoint[procI][pointI], nCurPoints))
            {
                // Add the point
                zonePoints[nCurPoints] = curPoint;

                // Record point mapping
                pointMap[pointI] = nCurPoints;
                nCurPoints++;
            }
            else
            {
                // This point has already been added so we will not add it again
                nDuplicatePoints++;

                // Lookup previous instance of the point
                const label pID = procPointMap[ownerPoint[procI][pointI]];

                // Set the map to point to the previous instance of the point
                pointMap[pointI] = pID;
            }
        }

        // Add faces from all processors
        const faceList& curProcFaces = procPatchFaces[procI];

        // Label face map for the current processor
        labelList faceMap(curProcFaces.size());

        forAll (curProcFaces, faceI)
        {
            // Renumber face into new points
            face curFace = curProcFaces[faceI];

            forAll (curFace, fI)
            {
                curFace[fI] = pointMap[curFace[fI]];
            }

            // Record the face into zone
            zoneFaces[nCurFaces] = curFace;
            faceMap[faceI] = nCurFaces;
            nCurFaces++;
        }

        if (procI == Pstream::myProcNo())
        {
            // Store point addressing
            pointToGlobalAddrPtr_.set(new labelList(pointMap));

            // Store face addressing
            faceToGlobalAddrPtr_.set(new labelList(faceMap));
        }
    }

    // Resize the points list
    zonePoints.resize(nCurPoints);

    // All points and faces are collected.  Make a patch
    globalPatchPtr_.set(new standAlonePatch(zoneFaces, zonePoints));

    if (debug)
    {
        InfoInFunction
            << "Finished calculating primitive patch" << nl
            << "    nDuplicatePoints: " << nDuplicatePoints << endl;
    }
}


void Foam::globalPolyPatch::calcGlobalMasterToCurrentProcPointAddr() const
{
    if (globalMasterToCurrentProcPointAddrPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set"
            << abort(FatalError);
    }

    globalMasterToCurrentProcPointAddrPtr_.set
    (
        new labelList(globalPatch().nPoints(), -1)
    );
    labelList& curMap = globalMasterToCurrentProcPointAddrPtr_();

    vectorField fzGlobalPoints(globalPatch().localPoints());

    // Set all slave points to zero because only the master order is used
    if (!Pstream::master())
    {
        fzGlobalPoints = vector::zero;
    }

    // Pass points to all procs
    reduce(fzGlobalPoints, sumOp<List<vector>>());

    // Now every proc has the master's list of FZ points
    // every proc must now find the mapping from their local FZ points to
    // the master FZ points

    const vectorField& fzLocalPoints = globalPatch().localPoints();

    const edgeList& fzLocalEdges = globalPatch().edges();

    const labelListList& fzPointEdges = globalPatch().pointEdges();

    scalarField minEdgeLength(fzLocalPoints.size(), GREAT);

    forAll(minEdgeLength, pI)
    {
        const labelList& curPointEdges = fzPointEdges[pI];

        forAll(curPointEdges, eI)
        {
            scalar Le = fzLocalEdges[curPointEdges[eI]].mag(fzLocalPoints);

            if (Le < minEdgeLength[pI])
            {
                minEdgeLength[pI] = Le;
            }
        }
    }

    forAll(fzGlobalPoints, globalPointI)
    {
        boolList visited(fzLocalPoints.size(), false);

        forAll(fzLocalPoints, procPointI)
        {
            if (!visited[procPointI])
            {
                visited[procPointI] = true;

                label nextPoint = procPointI;

                scalar curDist =
                    mag
                    (
                        fzLocalPoints[nextPoint]
                      - fzGlobalPoints[globalPointI]
                    );

                if (curDist < 1e-4*minEdgeLength[nextPoint])
                {
                    curMap[globalPointI] = nextPoint;
                    break;
                }

                label found = false;

                while (nextPoint != -1)
                {
                    const labelList& nextPointEdges =
                        fzPointEdges[nextPoint];

                    scalar minDist = GREAT;
                    label index = -1;
                    forAll(nextPointEdges, edgeI)
                    {
                        label curNgbPoint =
                            fzLocalEdges
                            [
                                nextPointEdges[edgeI]
                            ].otherVertex(nextPoint);

                        if (!visited[curNgbPoint])
                        {
                            visited[curNgbPoint] = true;

                            scalar curDist =
                                mag
                                (
                                    fzLocalPoints[curNgbPoint]
                                  - fzGlobalPoints[globalPointI]
                                );

                            if (curDist < 1e-4*minEdgeLength[curNgbPoint])
                            {
                                curMap[globalPointI] = curNgbPoint;
                                found = true;
                                break;
                            }
                            else if (curDist < minDist)
                            {
                                minDist = curDist;
                                index = curNgbPoint;
                            }
                        }
                    }

                    nextPoint = index;
                }

                if (found)
                {
                    break;
                }
            }
        }
    }

    forAll(curMap, globalPointI)
    {
        if (curMap[globalPointI] == -1)
        {
            FatalErrorInFunction
                << "point map is not correct!"
                << abort(FatalError);
        }
    }
}


void Foam::globalPolyPatch::calcInterp() const
{
    if (debug)
    {
        InfoInFunction
            << "Calculating patch interpolator"
            << endl;
    }

    if (interpPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set"
            << abort(FatalError);
    }

    interpPtr_.reset
    (
        new PrimitivePatchInterpolation<standAlonePatch>(globalPatch())
    );
}


void Foam::globalPolyPatch::calcLocalInterp() const
{
    if (debug)
    {
        InfoInFunction
            << "Calculating local patch interpolator"
            << endl;
    }

    if (localInterpPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set"
            << abort(FatalError);
    }

    localInterpPtr_.reset
    (
        new primitivePatchInterpolation(patch())
    );
}


void Foam::globalPolyPatch::check() const
{
    label patchIndex = mesh_.boundaryMesh().findPatchID(patchName_);

    if (patchIndex < 0)
    {
        FatalErrorInFunction
            << "Patch " << patchName_ << " not found."
            << abort(FatalError);
    }
}


void Foam::globalPolyPatch::clearOut() const
{
    globalPatchPtr_.clear();
    pointToGlobalAddrPtr_.clear();
    faceToGlobalAddrPtr_.clear();
    globalMasterToCurrentProcPointAddrPtr_.clear();
    interpPtr_.clear();
    localInterpPtr_.clear();
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
    patch_(mesh_.boundaryMesh()[mesh_.boundaryMesh().findPatchID(patchName_)]),
    displacementField_
    (
        dict.lookupOrDefault<word>("displacementField", "none")
    ),
    globalPatchPtr_(NULL),
    pointToGlobalAddrPtr_(NULL),
    faceToGlobalAddrPtr_(NULL),
    globalMasterToCurrentProcPointAddrPtr_(NULL),
    interpPtr_(NULL),
    localInterpPtr_(NULL)
{
    check();
}


Foam::globalPolyPatch::globalPolyPatch
(
    const polyPatch& patch,
    const word& displacementField
)
:
    mesh_(patch.boundaryMesh().mesh()),
    patchName_(patch.name()),
    patch_(mesh_.boundaryMesh()[mesh_.boundaryMesh().findPatchID(patchName_)]),
    displacementField_(displacementField),
    globalPatchPtr_(NULL),
    pointToGlobalAddrPtr_(NULL),
    faceToGlobalAddrPtr_(NULL),
    globalMasterToCurrentProcPointAddrPtr_(NULL),
    interpPtr_(NULL),
    localInterpPtr_(NULL)
{
    check();
}


Foam::autoPtr<Foam::globalPolyPatch> Foam::globalPolyPatch::New
(
    const dictionary& dict,
    const polyPatch& patch
)
{
    if (!patch.coupled())
    {
        return autoPtr<globalPolyPatch>
        (
            new globalPolyPatch(dict, patch)
        );
    }
    return autoPtr<globalPolyPatch>
    (
        new coupledGlobalPolyPatch(dict.subDict(patch.name()), patch)
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


const Foam::standAlonePatch& Foam::globalPolyPatch::globalPatch() const
{
    if (!globalPatchPtr_.valid())
    {
        calcGlobalPatch();

        if (displacementField_ != "none")
        {
            if (debug && mesh_.time().outputTime())
            {
                mkDir("VTK");
                globalPatchPtr_->writeVTK
                (
                    "VTK/"
                    + patch_.name() + '_'
                    + Foam::name(mesh_.time().timeIndex())
                );
            }
        }
    }

    return globalPatchPtr_();
}


const Foam::PrimitivePatchInterpolation<Foam::standAlonePatch>&
Foam::globalPolyPatch::interpolator() const
{
    if (!interpPtr_.valid())
    {
        calcInterp();
    }

    return interpPtr_();
}


const Foam::primitivePatchInterpolation&
Foam::globalPolyPatch::localInterpolator() const
{
    if (!localInterpPtr_.valid())
    {
        calcLocalInterp();
    }

    return localInterpPtr_();
}


const Foam::labelList& Foam::globalPolyPatch::pointToGlobalAddr() const
{
    if (!pointToGlobalAddrPtr_.valid())
    {
        calcGlobalPatch();
    }

    return pointToGlobalAddrPtr_();
}


const Foam::labelList& Foam::globalPolyPatch::faceToGlobalAddr() const
{
    if (!faceToGlobalAddrPtr_.valid())
    {
        calcGlobalPatch();
    }

    return faceToGlobalAddrPtr_();
}


const Foam::labelList&
Foam::globalPolyPatch::globalMasterToCurrentProcPointAddr() const
{
    if (!globalMasterToCurrentProcPointAddrPtr_.valid())
    {
        calcGlobalMasterToCurrentProcPointAddr();
    }

    return globalMasterToCurrentProcPointAddrPtr_();
}


void Foam::globalPolyPatch::updateMesh()
{
    clearOut();
}


void Foam::globalPolyPatch::movePoints()
{
    if (globalPatchPtr_.valid())
    {
        globalPatchPtr_->movePoints(patch_.points());
    }
}


void Foam::globalPolyPatch::movePoints(const pointField& pts)
{
    if (globalPatchPtr_.valid())
    {
        globalPatchPtr_->movePoints(pts);
    }
}


// ************************************************************************* //
