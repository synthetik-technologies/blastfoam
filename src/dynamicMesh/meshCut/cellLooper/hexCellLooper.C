/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "hexCellLooper.H"
#include "cellFeatures.H"
#include "polyMesh.H"
#include "cellModeller.H"
#include "plane.H"
#include "ListOps.H"
#include "meshTools.H"
#include "OFstream.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(hexCellLooper, 0);
addToRunTimeSelectionTable(cellLooper, hexCellLooper, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Starting from cut edge start walking.
bool Foam::hexCellLooper::walkHex
(
    const label celli,
    const label startFacei,
    const label startEdgeI,

    labelList& loop,
    scalarField& loopWeights
) const
{
    label facei = startFacei;

    label edgeI = startEdgeI;

    label cutI = 0;

    do
    {
        if (debug & 2)
        {
            Pout<< "    walkHex : inserting cut onto edge:" << edgeI
                << " vertices:" << mesh().edges()[edgeI] << endl;
        }

        // Store cut through edge. For now cut edges halfway.
        loop[cutI] = edgeToEVert(edgeI);
        loopWeights[cutI] = 0.5;
        cutI++;

        facei = meshTools::otherFace(mesh(), celli, facei, edgeI);

        const edge& e = mesh().edges()[edgeI];

        // Walk two edges further
        edgeI = meshTools::walkFace(mesh(), facei, edgeI, e.end(), 2);

        if (edgeI == startEdgeI)
        {
            break;
        }
    }
    while (true);

    // Checks.
    if (cutI > 4)
    {
        Pout<< "hexCellLooper::walkHex" << "Problem : cell:" << celli
            << " collected loop:";
        writeCuts(Pout, loop, loopWeights);
        Pout<< "loopWeights:" << loopWeights << endl;

        return false;
    }
    else
    {
        return true;
    }
}



void Foam::hexCellLooper::makeFace
(
    const labelList& loop,
    const scalarField& loopWeights,

    labelList& faceVerts,
    pointField& facePoints
) const
{
    facePoints.setSize(loop.size());
    faceVerts.setSize(loop.size());

    forAll(loop, cutI)
    {
        label cut = loop[cutI];

        if (isEdge(cut))
        {
            label edgeI = getEdge(cut);

            const edge& e = mesh().edges()[edgeI];

            const point& v0 = mesh().points()[e.start()];
            const point& v1 = mesh().points()[e.end()];

            facePoints[cutI] =
                loopWeights[cutI]*v1 + (1.0-loopWeights[cutI])*v0;
        }
        else
        {
            label vertI = getVertex(cut);

            facePoints[cutI] = mesh().points()[vertI];
        }

        faceVerts[cutI] = cutI;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::hexCellLooper::hexCellLooper(const polyMesh& mesh)
:
    geomCellLooper(mesh),
    hex_(*(cellModeller::lookup("hex")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::hexCellLooper::~hexCellLooper()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::hexCellLooper::cut
(
    const vector& refDir,
    const label celli,
    const boolList& vertIsCut,
    const boolList& edgeIsCut,
    const scalarField& edgeWeight,

    labelList& loop,
    scalarField& loopWeights
) const
{
    bool success = false;

    if (mesh().cellShapes()[celli].model() == hex_)
    {
        // Get starting edge. Note: should be compatible with way refDir is
        // determined.
        label edgeI = meshTools::cutDirToEdge(mesh(), celli, refDir);

        // Get any face using edge
        label face0;
        label face1;
        meshTools::getEdgeFaces(mesh(), celli, edgeI, face0, face1);

        // Walk circumference of hex, cutting edges only
        loop.setSize(4);
        loopWeights.setSize(4);

        success = walkHex(celli, face0, edgeI, loop, loopWeights);
    }
    else
    {
        success = geomCellLooper::cut
        (
            refDir,
            celli,
            vertIsCut,
            edgeIsCut,
            edgeWeight,

            loop,
            loopWeights
        );
    }

    if (debug)
    {
        if (loop.empty())
        {
            WarningInFunction
                << "could not cut cell " << celli << endl;

            fileName cutsFile("hexCellLooper_" + name(celli) + ".obj");

            Pout<< "hexCellLooper : writing cell to " << cutsFile << endl;

            OFstream cutsStream(cutsFile);

            meshTools::writeOBJ
            (
                cutsStream,
                mesh().cells(),
                mesh().faces(),
                mesh().points(),
                labelList(1, celli)
            );

            return false;
        }

        // Check for duplicate cuts.
        labelHashSet loopSet(loop.size());

        forAll(loop, elemI)
        {
            label elem = loop[elemI];

            if (loopSet.found(elem))
            {
                FatalErrorInFunction
                    << abort(FatalError);
            }
            loopSet.insert(elem);
        }


        face faceVerts(loop.size());
        pointField facePoints(loop.size());

        makeFace(loop, loopWeights, faceVerts, facePoints);

        if ((faceVerts.mag(facePoints) < small) || (loop.size() < 3))
        {
            FatalErrorInFunction
                << " on points:" << facePoints << endl
                << UIndirectList<point>(facePoints, faceVerts)()
                << abort(FatalError);
        }
    }
    return success;
}


// No shortcuts for cutting with arbitrary plane.
bool Foam::hexCellLooper::cut
(
    const plane& cutPlane,
    const label celli,
    const boolList& vertIsCut,
    const boolList& edgeIsCut,
    const scalarField& edgeWeight,

    labelList& loop,
    scalarField& loopWeights
) const
{
    return
        geomCellLooper::cut
        (
            cutPlane,
            celli,
            vertIsCut,
            edgeIsCut,
            edgeWeight,

            loop,
            loopWeights
        );
}


// ************************************************************************* //
