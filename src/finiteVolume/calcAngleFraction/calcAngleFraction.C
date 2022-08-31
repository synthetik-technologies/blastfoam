/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "calcAngleFraction.H"
#include "wedgePolyPatch.H"
#include "constants.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


Foam::scalar Foam::calcAngleFraction(const polyMesh& mesh)
{
    // Return 1 for non-wedge cases
    if (mesh.geometricD() == mesh.solutionD())
    {
        return 1.0;
    }

    scalar scale(1.0);
    List<vector> axes;
    scalar pi = Foam::constant::mathematical::pi;
    const polyBoundaryMesh& pbMesh(mesh.boundaryMesh());
    forAll(pbMesh, patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];
        if (isA<wedgePolyPatch>(patch))
        {
            const wedgePolyPatch& wedge
            (
                dynamicCast<const wedgePolyPatch>(patch)
            );
            bool found = false;

            // Only add one face per pair of wedges
            vector axis = cmptMag(wedge.axis());
            forAll(axes, i)
            {
                if (mag(axis - axes[i]) < 1e-6)
                {
                    found = true;
                }
            }

            // If the wedge pair has not yet been added, add it and
            // modify the scale factor
            // The range corresponds to the maximum angle for each
            // rotations (ie. 2*pi, pi, 0)
            if (!found)
            {
                axes.append(axis);
                scale *= acos(wedge.cosAngle());
            }
        }
    }
    // Integration over [0:pi, 0:p/2] results in
    // pi*(1 - cos(pi/2)) = pi
    return scale/pi;
}


// ************************************************************************* //
