if (adjustTimeStep)
{
    const scalar maxMDotCo
    (
        runTime.controlDict().lookupOrDefault("maxMDotCo", 0.1)
    );

    scalar deltaT = 1.2*runTime.deltaTValue();

    if (CoNum > small)
    {
        deltaT = min(deltaT, runTime.deltaTValue()*maxCo/CoNum);
    }

    if (mDotCoNum > small)
    {
        deltaT = min(maxMDotCo/mDotCoNum*runTime.deltaTValue(), deltaT);
    }

    deltaT = min(min(deltaT, fvModels.maxDeltaT()), maxDeltaT);

    runTime.setDeltaT(deltaT);

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}
