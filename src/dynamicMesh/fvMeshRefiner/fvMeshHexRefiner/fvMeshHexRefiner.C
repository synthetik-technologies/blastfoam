/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-05-2020 Synthetik Applied Technologies: |    Modified original
                            dynamicRefineBalanceBlastFvMesh class
                            to be more appilcable to compressible flows.
                            Improved compatibility with snappyHexMesh.
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

#include "fvMeshHexRefiner.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceInterpolate.H"
#include "volFields.H"
#include "polyTopoChange.H"
#include "surfaceFields.H"
#include "syncTools.H"
#include "pointFields.H"
#include "fvCFD.H"
#include "volPointInterpolation.H"
#include "pointMesh.H"
#include "cellSet.H"
#include "wedgePolyPatch.H"
#include "hexRef3D.H"
#include "parcelCloud.H"
#include "hexRefRefinementHistoryConstraint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMeshHexRefiner, 0);
    addToRunTimeSelectionTable(fvMeshRefiner, fvMeshHexRefiner, fvMesh);
    addToRunTimeSelectionTable(fvMeshRefiner, fvMeshHexRefiner, dictionary);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshHexRefiner::fvMeshHexRefiner(fvMesh& mesh)
:
    fvMeshRefiner(mesh)
{
    // Added refinement history decomposition constraint to keep all
    // cells with the same parent together
    {
        dictionary refinementHistoryDict("refinementHistory");
        refinementHistoryDict.add
        (
            "type",
            hexRefRefinementHistoryConstraint::typeName
        );
        balancer_.addConstraint("refinementHistory", refinementHistoryDict);
    }
}


Foam::fvMeshHexRefiner::fvMeshHexRefiner
(
    fvMesh& mesh,
    const dictionary& dict,
    const bool force,
    const bool read
)
:
    fvMeshRefiner(mesh, dict, force, read)
{
    // Added refinement history decomposition constraint to keep all
    // cells with the same parent together
    {
        dictionary refinementHistoryDict;
        refinementHistoryDict.add
        (
            "type",
            hexRefRefinementHistoryConstraint::typeName
        );
        balancer_.addConstraint("refinementHistory", refinementHistoryDict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshHexRefiner::~fvMeshHexRefiner()
{}

// ************************************************************************* //
