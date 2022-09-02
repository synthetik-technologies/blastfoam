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

autoPtr<timeIntegrator> integrator(timeIntegrator::New(mesh));
reactingCompressibleSystem fluid(mesh);
integrator->addSystem(fluid);

const volScalarField& p = fluid.p();
const volScalarField& T = fluid.T();
fluid.update();
