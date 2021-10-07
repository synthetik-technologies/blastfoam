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
#include "polyPatchID.H"

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
        InfoIn("void globalPolyPatch::calcGlobalPatch() const")
            << "Calculating primitive patch"
            << endl;
    }

    if (globalPatchPtr_ || pointToGlobalAddrPtr_ || faceToGlobalAddrPtr_)
    {
        FatalErrorIn
        (
            "void globalPolyPatch::calcGlobalPatch() const"
        )   << "primitive face zone patch and addressing already calculated"
            << abort(FatalError);
    }

    // Get patch
    polyPatchID patchID
    (
        patchName_,
        mesh_.boundaryMesh()
    );

    if (!patchID.active())
    {
        FatalErrorIn("void globalPolyPatch::calcGlobalPatch() const")
            << "Cannot find patch " << patchName_
            << abort(FatalError);
    }

    // Collect points and faces from all processors
    typedef List<point> pointList;
    typedef List<pointList> pointListList;

    pointListList procPatchPoints(Pstream::nProcs());
    faceListList procPatchFaces(Pstream::nProcs());

    // Add points and faces if the patch is not empty
    if (!mesh_.boundaryMesh()[patchID.index()].empty())
    {
        // Insert my points
        procPatchPoints[Pstream::myProcNo()] =
            mesh_.boundaryMesh()[patchID.index()].localPoints();

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
        FatalErrorIn("void globalPolyPatch::calcGlobalPatch() const")
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
    HashTable<label, point, Hash<point> > zonePointsSet(nZonePoints);
    label nDuplicatePoints = 0;

    label nCurPoints = 0;
    label nCurFaces = 0;

    // Collect all points and faces
    forAll (procPatchPoints, procI)
    {
        // Add points from all processors except self
        const pointList& curProcPoints = procPatchPoints[procI];

        // Label point map for the current processor
        labelList pointMap(curProcPoints.size());

        // Add points from all processors
        forAll (curProcPoints, pointI)
        {
            // Note: possible removal of duplicate points here
            // HJ, 28/Dec/2016
            // PC, 15/Dec/2017: the edge loops are wrong unless the duplicates
            // are removed, because the internal patch faces will not be
            // connected across the processor boundaries. So I will keep track
            // of the duplicates with a HashTable and remove them

            // Current point
            const point& curPoint = curProcPoints[pointI];

            // Check if the point has already been added
            HashTable<label, point, Hash<point> >::iterator iter =
                zonePointsSet.find(curPoint);

            if (iter == zonePointsSet.end())
            {
                // This is a new point; add it and record the index
                zonePointsSet.insert(curPoint, nCurPoints);

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
                const label pID = iter();

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
            pointToGlobalAddrPtr_ = new labelList(pointMap);

            // Store face addressing
            faceToGlobalAddrPtr_ = new labelList(faceMap);
        }
    }

    // Resize the points list
    zonePoints.resize(nCurPoints);

    // All points and faces are collected.  Make a patch
    globalPatchPtr_ = new standAlonePatch(zoneFaces, zonePoints);

    if (debug)
    {
        Info<< "    nDuplicatePoints: " << nDuplicatePoints << endl;
    }

    if (debug)
    {
        Info<< "void globalPolyPatch::calcGlobalPatch() const : "
            << "Finished calculating primitive patch"
            << endl;
    }
}


void Foam::globalPolyPatch::calcGlobalMasterToCurrentProcPointAddr() const
{
    if (globalMasterToCurrentProcPointAddrPtr_)
    {
        FatalErrorIn
        (
            "void globalPolyPatch::calcGlobalMasterToCurrentProcPointAddr() "
            "const"
        )   << "pointer already set"
            << abort(FatalError);
    }

    globalMasterToCurrentProcPointAddrPtr_ =
        new labelList(globalPatch().nPoints(), -1);
    labelList& curMap = *globalMasterToCurrentProcPointAddrPtr_;

    vectorField fzGlobalPoints(globalPatch().localPoints());

    // Pass points to all procs
    Pstream::scatterList<vector>(fzGlobalPoints);

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
            FatalErrorIn
            (
                type() + "::calcGlobalMasterToCurrentProcPointAddr()"
            )   << "point map is not correct!"
                << abort(FatalError);
        }
    }
}


void Foam::globalPolyPatch::calcInterp() const
{
    if (debug)
    {
        InfoIn("void globalPolyPatch::calcInterp() const")
            << "Calculating patch interpolator"
            << endl;
    }

    if (interpPtr_)
    {
        FatalErrorIn
        (
            "void globalPolyPatch::calcInterp() const"
        )   << "pointer already set"
            << abort(FatalError);
    }

    interpPtr_ =
        new PrimitivePatchInterpolation<standAlonePatch>(globalPatch());
}


void Foam::globalPolyPatch::check() const
{
    label patchIndex = mesh_.boundaryMesh().findPatchID(patchName_);

    if (patchIndex < 0)
    {
        FatalErrorIn("void globalPolyPatch::check() const")
            << "Patch " << patchName_ << " not found."
            << abort(FatalError);
    }
}


void Foam::globalPolyPatch::clearOut() const
{
    deleteDemandDrivenData(globalPatchPtr_);
    deleteDemandDrivenData(pointToGlobalAddrPtr_);
    deleteDemandDrivenData(faceToGlobalAddrPtr_);
    deleteDemandDrivenData(globalMasterToCurrentProcPointAddrPtr_);
    deleteDemandDrivenData(interpPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::globalPolyPatch::globalPolyPatch
(
    const word& patchName,
    const polyMesh& mesh
)
:
    mesh_(mesh),
    patchName_(patchName),
    patch_(mesh_.boundaryMesh()[mesh_.boundaryMesh().findPatchID(patchName_)]),
    globalPatchPtr_(NULL),
    pointToGlobalAddrPtr_(NULL),
    faceToGlobalAddrPtr_(NULL),
    globalMasterToCurrentProcPointAddrPtr_(NULL),
    interpPtr_(NULL)
{
    check();
}


// Construct from dictionary
Foam::globalPolyPatch::globalPolyPatch
(
    const dictionary& dict,
    const polyMesh& mesh
)
:
    mesh_(mesh),
    patchName_(dict.lookup("patch")),
    patch_(mesh_.boundaryMesh()[mesh_.boundaryMesh().findPatchID(patchName_)]),
    globalPatchPtr_(NULL),
    pointToGlobalAddrPtr_(NULL),
    faceToGlobalAddrPtr_(NULL),
    globalMasterToCurrentProcPointAddrPtr_(NULL),
    interpPtr_(NULL)
{
    check();
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
    if (!globalPatchPtr_)
    {
        calcGlobalPatch();
    }

    return *globalPatchPtr_;
}


const Foam::PrimitivePatchInterpolation<Foam::standAlonePatch>&
Foam::globalPolyPatch::interpolator() const
{
    if (!interpPtr_)
    {
        calcInterp();
    }

    return *interpPtr_;
}


const Foam::labelList& Foam::globalPolyPatch::pointToGlobalAddr() const
{
    if (!pointToGlobalAddrPtr_)
    {
        calcGlobalPatch();
    }

    return *pointToGlobalAddrPtr_;
}


const Foam::labelList& Foam::globalPolyPatch::faceToGlobalAddr() const
{
    if (!faceToGlobalAddrPtr_)
    {
        calcGlobalPatch();
    }

    return *faceToGlobalAddrPtr_;
}


const Foam::labelList&
Foam::globalPolyPatch::globalMasterToCurrentProcPointAddr() const
{
    if (!globalMasterToCurrentProcPointAddrPtr_)
    {
        calcGlobalMasterToCurrentProcPointAddr();
    }

    return *globalMasterToCurrentProcPointAddrPtr_;
}


void Foam::globalPolyPatch::updateMesh()
{
    clearOut();
}


void Foam::globalPolyPatch::movePoints(const pointField& p)
{
    if (globalPatchPtr_)
    {
        globalPatchPtr_->movePoints(p);
    }
}


// ************************************************************************* //
