//- Set the displacement references for mapping of moving patches
forAll(fluidRegions, i)
{
    setMappedPatchDisplacement(fluidRegions[i], "pointDisplacement");
}
