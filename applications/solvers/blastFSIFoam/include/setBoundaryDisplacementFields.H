//- Set the displacement references for mapping of moving patches
forAll(fluidRegions, i)
{
    mappedPatchSelector::replaceMappedWithMappedMoving(fluidRegions[i]);
    mappedPointPatchSelector::replaceMappedWithMappedMoving
    (
        fluidRegions[i]
    );

    forAll(solidRegions, j)
    {
        mappedPatchSelector::setMappedPatchDisplacement
        (
            fluidRegions[i],
            solidRegions[j]
        );
        mappedPointPatchSelector::setMappedPatchDisplacement
        (
            fluidRegions[i],
            solidRegions[j]
        );
    }
}
