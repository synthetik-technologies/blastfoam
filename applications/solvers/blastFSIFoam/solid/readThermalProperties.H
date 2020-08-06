thermalStresses[i] =
    thermophysicalProperties[i].lookupType<Switch>("thermalStress");

if (thermalStresses[i])
{
    IOobject alphaIO
    (
        "alpha",
        runTime.timeName(0),
        solidRegions[i],
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    const dictionary& alphaDict(thermophysicalProperties[i].subDict("alpha"));
    word alphaType(alphaDict.lookup("type"));
    if (alphaType == "uniform")
    {
        scalar alphaValue(readScalar(alphaDict.lookup("value")));

        alphas.set
        (
            i,
            new volScalarField
            (
                alphaIO,
                solidRegions[i],
                dimensionedScalar
                (
                    dimensionSet(0, 0, 0, -1, 0),
                    alphaValue
                )
            )
        );

    }
    else if (alphaType == "field")
    {
        alphaIO.readOpt() = IOobject::MUST_READ;

        alphas.set
        (
            i,
            new volScalarField
            (
                alphaIO,
                solidRegions[i]
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Valid type entries are uniform or field for alpha"
            << abort(FatalError);
    }
}
