    scalar CoNum = -great;

    forAll(fluidRegions, regionI)
    {
        CoNum = max
        (
            fluids[regionI].CoNum(),
            CoNum
        );
    }
