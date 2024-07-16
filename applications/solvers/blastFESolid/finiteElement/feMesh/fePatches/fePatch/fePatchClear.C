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

#include "fePatch.H"
#include "demandDrivenData.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fePatch::clearGeom()
{
    if (debug)
    {
        InfoInFunction << "Clearing geometric data" << endl;
    }

    deleteDemandDrivenData(shapesPtr_);
    deleteDemandDrivenData(dshapesPtr_);
    deleteDemandDrivenData(Js_);
    deleteDemandDrivenData(invJs_);
    deleteDemandDrivenData(Ws_);
    deleteDemandDrivenData(Bs_);
    deleteDemandDrivenData(nodeNormalsPtr_);
}


void Foam::fePatch::clearTopology()
{
    if (debug)
    {
        InfoInFunction << "Clearing patch addressing" << endl;
    }

    deleteDemandDrivenData(nodeEdgesPtr_);
    deleteDemandDrivenData(nodeFacesPtr_);
}


void Foam::fePatch::clearPatchMeshAddr()
{
    if (debug)
    {
        InfoInFunction << "Clearing patch-mesh addressing" << endl;
    }

    deleteDemandDrivenData(meshNodesPtr_);
    deleteDemandDrivenData(meshNodeMapPtr_);
}


void Foam::fePatch::clearOut()
{
    clearGeom();
    clearTopology();
    clearPatchMeshAddr();
}


// ************************************************************************* //
