mechanicalProperties.set
(
    i,
    new IOdictionary
    (
        IOobject
        (
            "mechanicalProperties",
            runTime.constant(),
            solidRegions[i],
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    )
);

const dictionary& EDict(mechanicalProperties[i].subDict("E"));
word EType(EDict.lookup("type"));

IOobject EHeader
(
    "E",
    runTime.timeName(0),
    solidRegions[i],
    IOobject::NO_READ,
    IOobject::NO_WRITE
);

if (EType == "uniform")
{
    scalar rhoEValue(readScalar(EDict.lookup("value")));

    Es.set
    (
        i,
        new volScalarField
        (
            EHeader,
            solidRegions[i],
            dimensionedScalar
            (
                dimPressure,
                rhoEValue
            )
        )
    );
}
else if (EType == "field")
{
    EHeader.readOpt() = IOobject::MUST_READ;

    Es.set
    (
        i,
        new volScalarField
        (
            EHeader,
            solidRegions[i]
        )
    );
}
else
{
    FatalErrorInFunction
        << "Valid type entries are uniform or field for E"
        << abort(FatalError);
}

IOobject nuIO
(
    "nu",
    runTime.timeName(0),
    solidRegions[i],
    IOobject::NO_READ,
    IOobject::NO_WRITE
);

const dictionary& nuDict(mechanicalProperties[i].subDict("nu"));
word nuType(nuDict.lookup("type"));

if (nuType == "uniform")
{
    scalar nuValue(readScalar(nuDict.lookup("value")));
    nus.set
    (
        i,
        new volScalarField
        (
            nuIO,
            solidRegions[i],
            dimensionedScalar
            (
                dimless,
                nuValue
            )
        )
    );
}
else if (nuType == "field")
{
    nuIO.readOpt() = IOobject::MUST_READ;
    nus.set
    (
        i,
        new volScalarField
        (
            nuIO,
            solidRegions[i]
        )
    );
}
else
{
    FatalErrorInFunction
        << "Valid type entries are uniform or field for nu"
        << abort(FatalError);
}

Info<< "Normalising E : E/rho\n" << endl;
Es[i] /= solidThermos[i].rho();

Info<< "Calculating Lame's coefficients\n" << endl;

const volScalarField& nu = nus[i];
const volScalarField& E = Es[i];

mus.set(i, new volScalarField(E/(2.0*(1.0 + nu))));
lambdas.set(i, new volScalarField(nu*E/((1.0 + nu)*(1.0 - 2.0*nu))));

planeStresses[i] = mechanicalProperties[i].lookupType<Switch>("planeStress");

if (planeStresses[i])
{
    Info<< "Plane Stress\n" << endl;
    lambdas[i] = nu*E/((1.0 + nu)*(1.0 - nu));
}
else
{
    Info<< "Plane Strain\n" << endl;
}
