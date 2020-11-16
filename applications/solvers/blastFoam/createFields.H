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
autoPtr<phaseCompressibleSystem> fluid
(
    phaseCompressibleSystem::New(mesh)
);
integrator->addSystem(fluid());

const volScalarField& p = fluid->p();
const volScalarField& T = fluid->T();
const surfaceScalarField& phi = fluid->phi();
fluid->update();

