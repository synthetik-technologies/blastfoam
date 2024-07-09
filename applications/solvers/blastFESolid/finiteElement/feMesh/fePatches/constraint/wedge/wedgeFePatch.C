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

#include "wedgeFePatch.H"
#include "pointConstraint.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wedgeFePatch, 0);

    // Add the patch constructor functions to the hash tables
    addToRunTimeSelectionTable
    (
        fePatch,
        wedgeFePatch,
        polyPatch
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wedgeFePatch::wedgeFePatch
(
    const polyPatch& patch,
    const feBoundaryMesh& bm
)
:
    fePatch(patch, bm),
    wedgePolyPatch_(refCast<const wedgePolyPatch>(patch))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wedgeFePatch::applyConstraint
(
    const label nodei,
    pointConstraint& pc
) const
{
    pc.applyConstraint(wedgePolyPatch_.n());
}


// ************************************************************************* //
