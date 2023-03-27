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

fvTimeIntegrator integrator(mesh);
reactingCompressibleSystem fluid(mesh);
integrator.addSystem(fluid);

const volScalarField& p = fluid.p();
const volScalarField& T = fluid.T();
fluid.update();
