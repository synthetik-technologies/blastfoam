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
    remapped[j] = fluidRegions[i].refine(false);
    needRemap = needRemap || remapped[j];
    j++;
}
forAll(solidRegions, i)
{
    // Already cleared mapped patches
    remapped[j] = solidRegions[i].refine(false);
    needRemap = needRemap || remapped[j];
    j++;
}
