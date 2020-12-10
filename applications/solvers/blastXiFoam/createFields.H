
// Read gravitational acceleration
Info<< "\nReading g" << endl;
uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    dimensionedVector(dimAcceleration, Zero)
);

autoPtr<timeIntegrator> integrator(timeIntegrator::New(mesh));
psiuCompressibleSystem fluid(mesh);
integrator->addSystem(fluid);

psiuReactionThermo& thermo = fluid.thermo();
basicCombustionMixture& composition = thermo.composition();
compressible::momentumTransportModel& turbulence = fluid.turbulence();
fluidThermophysicalTransportModel& thermophysicalTransport =
    fluid.thermophysicalTransport();

const volScalarField& rho = fluid.rho();
const volScalarField& p = thermo.p();
const volVectorField& U = fluid.U();
const surfaceScalarField& phi = fluid.rhoPhi();


volScalarField& b = composition.Y("b");
Info<< "min(b) = " << min(b).value() << endl;


Info<< "Creating field Xi\n" << endl;
volScalarField Xi
(
    IOobject
    (
        "Xi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Creating the unstrained laminar flame speed\n" << endl;
autoPtr<laminarFlameSpeed> unstrainedLaminarFlameSpeed
(
    laminarFlameSpeed::New(thermo)
);


Info<< "Reading strained laminar flame speed field Su\n" << endl;
volScalarField Su
(
    IOobject
    (
        "Su",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

dimensionedScalar SuMin = 0.01*Su.average();
dimensionedScalar SuMax = 4*Su.average();

Info<< "Calculating turbulent flame speed field St\n" << endl;
volScalarField St
(
    IOobject
    (
        "St",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Xi*Su
);
