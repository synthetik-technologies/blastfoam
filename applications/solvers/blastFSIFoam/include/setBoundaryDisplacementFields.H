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

Global
    setInitialMultiRegionDeltaT

Description
    Set the initial timestep for the CHT MultiRegion solver.

\*---------------------------------------------------------------------------*/


// Clear mapped boundaries if one region has been balanced
// Balanced meshes have already had their maps cleared
forAll(fluidRegions, i)
{
    if (fluidRegions[i].foundObject<volVectorField>("cellDisplacement"))
    {
        const volVectorField& cd =
            fluidRegions[i].lookupObject<volVectorField>("cellDisplacement");

        forAll(fluidRegions[i].boundaryMesh(), patchi)
        {
            if (isA<mappedMovingWallFvPatch>(fluidRegions[i].boundary()[patchi]))
            {
                polyBoundaryMesh& pbMesh =
                    const_cast<polyBoundaryMesh&>
                    (
                        fluidRegions[i].boundaryMesh()
                    );
                refCast<mappedMovingPatchBase>(pbMesh[patchi]).setOffsets(cd);
            }
        }
    }
}

// Solid meshes do not move so there is no displacement needed

// ************************************************************************* //
