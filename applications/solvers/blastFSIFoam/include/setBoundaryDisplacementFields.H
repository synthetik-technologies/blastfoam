//- Set the displacement references for mapping of moving patches
forAll(fluidRegions, i)
{
    forAll(solidRegions, j)
    {
        mappedPatchSelector::setMappedPatchDisplacement
        (
            fluidRegions[i],
            solidRegions[j],
            pointDs[j]
        );
        mappedPointPatchSelector::setMappedPatchDisplacement
        (
            fluidRegions[i],
            solidRegions[j],
            pointDs[j]
        );
    }
}
