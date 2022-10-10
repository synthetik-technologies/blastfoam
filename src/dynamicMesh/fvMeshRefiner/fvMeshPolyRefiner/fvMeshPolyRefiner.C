/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "fvMeshPolyRefiner.H"
#include "polyTopoChange.H"
#include "parcelCloud.H"
#include "prismatic2DRefinement.H"
#include "polyhedralRefinement.H"
#include "polyRefinementConstraint.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMeshPolyRefiner, 0);
    addToRunTimeSelectionTable(fvMeshRefiner, fvMeshPolyRefiner, fvMesh);
    addToRunTimeSelectionTable(fvMeshRefiner, fvMeshPolyRefiner, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshPolyRefiner::fvMeshPolyRefiner(fvMesh& mesh)
:
    fvMeshRefiner(mesh)
{
    // Added refinement history decomposition constraint to keep all
    // cells with the same parent together
    {
        dictionary refinementHistoryDict;
        refinementHistoryDict.add
        (
            "type",
            polyRefinementConstraint::typeName
        );
        balancer_.addConstraint("refinementHistory", refinementHistoryDict);
    }
}


Foam::fvMeshPolyRefiner::fvMeshPolyRefiner
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
        dictionary refinementHistoryDict("refinementHistory");
        refinementHistoryDict.add
        (
            "type",
            polyRefinementConstraint::typeName
        );
        balancer_.addConstraint("refinementHistory", refinementHistoryDict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshPolyRefiner::~fvMeshPolyRefiner()
{}

// ************************************************************************* //
