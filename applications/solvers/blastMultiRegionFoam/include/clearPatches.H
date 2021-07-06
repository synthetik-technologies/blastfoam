// Clear mapped boundaries if one region has been remapped
// Balanced meshes have already had their maps cleared
if (needRemap)
{
    forAll(fluidRegions, i)
    {
        mappedPatchSelector::clearMappedPatches(fluidRegions[i]);
        mappedPointPatchSelector::clearMappedPatches(fluidRegions[i]);
    }
    forAll(solidRegions, i)
    {
        mappedPatchSelector::clearMappedPatches(solidRegions[i]);
        mappedPointPatchSelector::clearMappedPatches(solidRegions[i]);
    }
}

// ************************************************************************* //
