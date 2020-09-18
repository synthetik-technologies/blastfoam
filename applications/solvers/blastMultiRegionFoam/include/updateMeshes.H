//- Update Meshes and check if mesh moved
j = 0;
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
