/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
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

#include "dynamicRefineBalancedBlastFvMesh.H"
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
#include "RefineBalanceMeshObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicRefineBalancedBlastFvMesh, 0);
    addToRunTimeSelectionTable
    (
        dynamicBlastFvMesh,
        dynamicRefineBalancedBlastFvMesh,
        IOobject
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicRefineBalancedBlastFvMesh::dynamicRefineBalancedBlastFvMesh
(
    const IOobject& io
)
:
    dynamicRefineBlastFvMesh(io),
    balancer_
    (
        *this,
        dynamicMeshDict().optionalSubDict("loadBalance")
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicRefineBalancedBlastFvMesh::~dynamicRefineBalancedBlastFvMesh()
{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicRefineBalancedBlastFvMesh::update()
{
    return false;
}


bool Foam::dynamicRefineBalancedBlastFvMesh::refine(const bool)
{
    bool hasChanged = dynamicRefineBlastFvMesh::refine();
    balancer_.read(dynamicMeshDict().optionalSubDict("loadBalance"));

    if (balancer_.canBalance() && hasChanged)
    {
        balancer_.distribute();
    }

    return hasChanged;
}


// ************************************************************************* //
