//- Update Meshes and check if mesh moved
forAll(fluidRegions, i)
{
    needRemap = needRemap || fluidRegions[i].update();
}
forAll(solidRegions, i)
{
    needRemap = needRemap || solidRegions[i].update();
}

#include "clearPatches.H"
