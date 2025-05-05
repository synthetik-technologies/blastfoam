// Read gravitational acceleration
Info<< "\nReading g" << endl;
uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    dimensionedVector(dimAcceleration, Zero)
);

IOdictionary physicalProperties
(
    IOobject
    (
        ::Foam::physicalProperties::typeName,
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);
fvTimeIntegrator integrator(mesh);
reactingCompressibleSystem fluid(physicalProperties, mesh);
integrator.addSystem(fluid);

const volScalarField& p = fluid.p();
const volScalarField& T = fluid.T();
fluid.update();

const ::Foam::fvModels& fvModels = integrator.models();
