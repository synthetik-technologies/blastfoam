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

Info<< "Reading phaseProperties dictionary\n" << endl;
IOdictionary phaseProperties
(
     IOobject
     (
          "phaseProperties",
          runTime.constant(),
          mesh,
          IOobject::MUST_READ,
          IOobject::NO_WRITE
     )
);

autoPtr<phaseCompressibleSystem> fluid
(
    phaseCompressibleSystem::New(mesh, phaseProperties)
);
autoPtr<timeIntegrator> integrator(timeIntegrator::New(mesh));
integrator->addSystem(fluid());

const volScalarField& p = fluid->p();
const volScalarField& T = fluid->T();
const surfaceScalarField& phi = fluid->phi();
fluid->update();

