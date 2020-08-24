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

//- Update internal error fields of all regions
forAll(fluidRegions, i)
{
    fluidRegions[i].updateError();
}
forAll(solidRegions, i)
{
    solidRegions[i].updateError();
}

//- Update error field boundaries
//  Done after all region errors have been updated to make sure
//  all fields are up to date
forAll(fluidRegions, i)
{
    fluidRegions[i].updateErrorBoundaries();
}
forAll(solidRegions, i)
{
    solidRegions[i].updateErrorBoundaries();
}

//- Update Meshes and check if balancing has occurred
boolList remapped(fluidRegions.size() + solidRegions.size(), false);
bool needRemap = false;
label j = 0;
forAll(fluidRegions, i)
{
    // Already cleared mapped patches
    remapped[j] = fluidRegions[i].update();
    needRemap = needRemap || remapped[j];
    j++;
}
forAll(solidRegions, i)
{
    // Already cleared mapped patches
    remapped[j] = solidRegions[i].update();
    needRemap = needRemap || remapped[j];
    j++;
}

// Clear mapped boundaries if one region has been remapped
// Balanced meshes have already had their maps cleared
if (needRemap)
{
    j = 0;
    forAll(fluidRegions, i)
    {
        if (!remapped[j])
        {
            forAll(fluidRegions[i].boundaryMesh(), patchi)
            {
                if (isA<mappedWallFvPatch>(fluidRegions[i].boundary()[patchi]))
                {
                    polyBoundaryMesh& pbMesh =
                        const_cast<polyBoundaryMesh&>
                        (
                            fluidRegions[i].boundaryMesh()
                        );
                    refCast<mappedPatchBase>(pbMesh[patchi]).clearOut();
                }
                if (isA<mappedMovingWallFvPatch>(fluidRegions[i].boundary()[patchi]))
                {
                    polyBoundaryMesh& pbMesh =
                        const_cast<polyBoundaryMesh&>
                        (
                            fluidRegions[i].boundaryMesh()
                        );
                    refCast<mappedMovingPatchBase>(pbMesh[patchi]).clearOut();
                }
            }
        }
        j++;
    }
    forAll(solidRegions, i)
    {
        if (!remapped[j])
        {
            forAll(solidRegions[i].boundaryMesh(), patchi)
            {
                if (isA<mappedWallFvPatch>(solidRegions[i].boundary()[patchi]))
                {
                    polyBoundaryMesh& pbMesh =
                        const_cast<polyBoundaryMesh&>
                        (
                            solidRegions[i].boundaryMesh()
                        );
                    refCast<mappedPatchBase>(pbMesh[patchi]).clearOut();
                }
                if (isA<mappedMovingWallFvPatch>(solidRegions[i].boundary()[patchi]))
                {
                    polyBoundaryMesh& pbMesh =
                        const_cast<polyBoundaryMesh&>
                        (
                            solidRegions[i].boundaryMesh()
                        );
                    refCast<mappedMovingPatchBase>(pbMesh[patchi]).clearOut();
                }
            }
        }
        j++;
    }
}

// ************************************************************************* //
