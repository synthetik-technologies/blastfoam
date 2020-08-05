thermalStresses[i] =
    thermophysicalProperties[i].lookupType<Switch>("thermalStress");

threeKalphas.set
(
    i,
    new volScalarField
    (
        IOobject
        (
            "threeKalpha",
            runTime.timeName(),
            solidRegions[i],
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        solidRegions[i],
        dimensionedScalar(dimensionSet(0, 2, -2 , -1, 0), 0)
    )
);


DTs.set
(
    i,
    new volScalarField
    (
        IOobject
        (
            "DT",
            runTime.timeName(),
            solidRegions[i],
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        solidRegions[i],
        dimensionedScalar(dimensionSet(0, 2, -1 , 0, 0), 0)
    )
);


if (thermalStresses[i])
{
    IOobject CIO
    (
        "C",
        runTime.timeName(0),
        solidRegions[i],
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    const dictionary& CDict(thermophysicalProperties[i].subDict("C"));
    word CType(CDict.lookup("type"));
    if (CType == "uniform")
    {
        scalar CValue(readScalar(CDict.lookup("value")));

        Cs.set
        (
            i,
            new volScalarField
            (
                CIO,
                solidRegions[i],
                dimensionedScalar
                (
                    dimensionSet(0, 2, -2 , -1, 0),
                    CValue
                )
            )
        );

    }
    else if (CType == "field")
    {
        CIO.readOpt() = IOobject::MUST_READ;

        Cs.set
        (
            i,
            new volScalarField
            (
                CIO,
                solidRegions[i]
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Valid type entries are uniform or field for C"
            << abort(FatalError);
    }

    IOobject rhoKIO
    (
        "k",
        runTime.timeName(0),
        solidRegions[i],
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    const dictionary& kDict(thermophysicalProperties[i].subDict("k"));
    word kType(kDict.lookup("type"));
    if (kType == "uniform")
    {
        scalar rhoKValue(readScalar(kDict.lookup("value")));

        ks.set
        (
            i,
            new volScalarField
            (
                rhoKIO,
                solidRegions[i],
                dimensionedScalar
                (
                    dimensionSet(1, 1, -3 , -1, 0),
                    rhoKValue
                )
            )
        );

    }
    else if (kType == "field")
    {
        rhoKIO.readOpt() = IOobject::MUST_READ;

        ks.set
        (
            i,
            new volScalarField
            (
                rhoKIO,
                solidRegions[i]
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Valid type entries are uniform or field for K"
            << abort(FatalError);
    }
    ks[i] /= solidThermos[i].rho();

    alphas.set
    (
        i,
        new volScalarField
        (
            solidThermos[i].kappa()/solidThermos[i].Cp()
        )
    );


    Info<< "Calculating thermal coefficients\n" << endl;

    threeKalphas[i] = threeKs[i]*alphas[i];
    DTs[i] = ks[i]/Cs[i];
}
