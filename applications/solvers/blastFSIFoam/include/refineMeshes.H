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
bool needRemap = false;
forAll(fluidRegions, i)
{
    // Already cleared mapped patches
    needRemap = needRemap || fluidRegions[i].refine(false);
}
forAll(solidRegions, i)
{
    // Already cleared mapped patches
    needRemap = needRemap || solidRegions[i].refine(false);
}

#include "clearPatches.H"
