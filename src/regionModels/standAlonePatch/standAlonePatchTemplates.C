/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "standAlonePatch.H"
#include "vtkWritePolyData.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Patch>
Foam::standAlonePatch Foam::standAlonePatch::createGlobalPatch
(
    const Patch& patch
)
{
    if (!Pstream::parRun())
    {
        return standAlonePatch(patch, patch.points());
    }

    List<List<point>> gPoints(Pstream::nProcs());
    faceListList gFaces(Pstream::nProcs());

    // Insert my points
    gPoints[Pstream::myProcNo()] = patch.points();

    // Insert my faces
    gFaces[Pstream::myProcNo()] = patch;

    // Communicate points
    Pstream::gatherList(gPoints);
    Pstream::scatterList(gPoints);

    // Communicate faces
    Pstream::gatherList(gFaces);
    Pstream::scatterList(gFaces);

    label nPoints = 0;
    label nFaces = 0;
    forAll(gPoints, proci)
    {
        nPoints += gPoints[proci].size();
        nFaces += gFaces[proci].size();
    }
    pointField ps(nPoints);
    faceList fs(nFaces);

    label pi = 0;
    label fi = 0;
    forAll(gPoints, proci)
    {
        const label start = pi;
        forAll(gPoints[proci], pj)
        {
            ps[pi++] = gPoints[proci][pj];
        }
        forAll(gFaces[proci], fj)
        {
            fs[fi] = gFaces[proci][fj];
            face& f = fs[fi];
            forAll(f, fpi)
            {
                f[fpi] = gFaces[proci][fj][fpi] + start;
            }
            fi++;
        }
    }
    return standAlonePatch(move(fs), move(ps));
}

// ************************************************************************* //
