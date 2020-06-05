    scalar CoNum = -great;

    forAll(fluidRegions, regionI)
    {
        CoNum = max
        (
            compressibleCourantNo
            (
                fluidRegions[regionI],
                runTime,
                fluids[regionI]
            ),
            CoNum
        );
    }
