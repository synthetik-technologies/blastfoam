    scalar CoNum = 0.0;

    forAll(fluidRegions, regionI)
    {
        CoNum = max
        (
            fluids[regionI].CoNum(),
            CoNum
        );
    }
