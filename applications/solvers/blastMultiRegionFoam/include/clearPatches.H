// Clear mapped boundaries if one region has been remapped
// Balanced meshes have already had their maps cleared
if (needRemap)
{
    j = 0;
    forAll(fluidRegions, i)
    {
        if (!remapped[j])
        {
            clearMappedPatches(fluidRegions[i]);
        }
        j++;
    }
    forAll(solidRegions, i)
    {
        if (!remapped[j])
        {
            clearMappedPatches(solidRegions[i]);
        }
        j++;
    }
}

// ************************************************************************* //
