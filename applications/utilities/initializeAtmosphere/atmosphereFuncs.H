namespace Foam
{

void tableAtmosphere
(
    const dictionary& dict,
    const uniformDimensionedVectorField& g,
    const vector& dir,
    const volScalarField& h,
    volScalarField& p,
    volScalarField& rho
)
{
    dictionary pTableDict(dict.subDict("pTable"));
    lookupTable1D pTable
    (
        pTableDict.lookup("file"),
        pTableDict.lookupOrDefault<word>("pMod", "none"),
        pTableDict.lookupOrDefault<word>("hMod", "none")
    );
    dictionary TTableDict(dict.subDict("TTable"));
    lookupTable1D TTable
    (
        pTableDict.lookup("file"),
        pTableDict.lookupOrDefault<word>("TMod", "none"),
        pTableDict.lookupOrDefault<word>("hMod", "none")
    );

    scalar R = 8314.32/28.97;
    forAll(rho, celli)
    {
        p[celli] = pTable.lookup(h[celli]);
        scalar T = TTable.lookup(h[celli]);
        rho[celli] = p[celli]/T/R;
    }

    forAll(rho.boundaryField(), patchi)
    {
        forAll(rho.boundaryField()[patchi], facei)
        {
            p.boundaryFieldRef()[patchi][facei] =
                pTable.lookup(h.boundaryField()[patchi][facei]);
            scalar T = TTable.lookup(h.boundaryField()[patchi][facei]);
            rho.boundaryFieldRef()[patchi][facei] =
                p.boundaryField()[patchi][facei]/T/R;
        }
    }
}


void standardAtmosphere76i
(
    const scalar& H,
    const scalar& p0,
    const scalar& rho0,
    const scalar& T0,
    scalar& p,
    scalar& rho
)
{
    //- Scale to km
    scalar h = H/1000.0;

    //- Gravitational acceleration * molecular weight / ideal gas constant
    scalar GMR = 9.80665*28.97/8.31432;

    scalar gradT, x;

    scalar T = T0;
    scalar Told = T0;
    p = 1.0;

    // 0 - 11 km
    if (h > 0)
    {
        x = min(h, 11.0);
        gradT = -6.5;
        T += gradT*x;
        p *= exp(log(T0/T)*GMR/gradT);
        h -= 11.0;
    }

    // 11 - 20 km
    if (h > 0)
    {
        x = min(h, 9.0);
        Told = T;
        T = Told;
        p *= exp(-GMR*x/T);
        h -= 9.0;
    }

    // 20 - 32 km
    if (h > 0)
    {
        x = min(h, 12.0);
        gradT = 1.0;
        Told = T;
        T += gradT*x;
        p *= exp(log(Told/T)*GMR/gradT);
        h -= 12.0;
    }

    // 32 - 47 km
    if (h > 0)
    {
        x = min(h, 15.0);
        gradT = 2.8;
        Told = T;
        T += gradT*x;
        p *= exp(log(Told/T)*GMR/gradT);
        h -= 15.0;
    }

    // 47 - 51 km
    if (h > 0)
    {
        x = min(h, 4.0);
        Told = T;
        p *= exp(-GMR*x/Told);
        h -= 4.0;
    }

    // 51 - 71 km
    if (h > 0)
    {
        x = min(h, 20.0);
        gradT = -2.8;
        Told = T;
        T += gradT*x;
        p *= exp(log(Told/T)*GMR/gradT);
        h -= 20.0;
    }

    // 71 - 84.852 km
    if (h > 0)
    {
        x = min(h, 13.852);
        gradT = -2.0;
        Told = T;
        T += gradT*x;
        p *= exp(log(Told/T)*GMR/gradT);
    }

    rho = rho0*p*T0/T;
    p *= p0;
}


void standardAtmosphere76
(
    const dictionary& dict,
    const uniformDimensionedVectorField& g,
    const vector& dir,
    const volScalarField& h,
    volScalarField& p,
    volScalarField& rho
)
{
    scalar p0 = 101325.0;
    scalar rho0 = 1.225;
    scalar T0 = 288.15;
    forAll(rho, celli)
    {
        standardAtmosphere76i
        (
            h[celli],
            p0,
            rho0,
            T0,
            p[celli],
            rho[celli]
        );
    }

    forAll(rho.boundaryField(), patchi)
    {
        forAll(rho.boundaryField()[patchi], facei)
        {
            standardAtmosphere76i
            (
                h.boundaryField()[patchi][facei],
                p0,
                rho0,
                T0,
                p.boundaryFieldRef()[patchi][facei],
                rho.boundaryFieldRef()[patchi][facei]
            );
        }
    }
}


void simpleAtmosphere
(
    const dictionary& dict,
    const uniformDimensionedVectorField& g,
    const vector& dir,
    const volScalarField& h,
    volScalarField& p,
    volScalarField& rho
)
{
    dimensionedScalar pRef("pRef", dimPressure, dict);
    p = pRef + (g & dir)*rho*h;
    p.correctBoundaryConditions();
}

}
