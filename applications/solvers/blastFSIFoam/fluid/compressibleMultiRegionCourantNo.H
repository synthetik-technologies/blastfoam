    scalar CoNum = -great;

    forAll(fluidRegions, regionI)
    {
        scalar regionCoNum = fluids[regionI].CoNum();

        Info<< fluidRegions[regionI].name() << " max Courant Number = "
            << regionCoNum << endl;
        CoNum = max
        (
            regionCoNum,
            CoNum
        );
    }
