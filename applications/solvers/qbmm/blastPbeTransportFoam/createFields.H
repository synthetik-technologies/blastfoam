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

volScalarField& p = fluid->p();
const volVectorField& U = fluid->U();
const surfaceScalarField& phi = fluid->phi();
fluid->update();


Info<< "Reading populationBalanceProperties\n" << endl;

IOdictionary populationBalanceProperties
(
    IOobject
    (
        "populationBalanceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


autoPtr<ODEPopulationBalanceModel> populationBalance
(
    ODEPopulationBalanceModel::New
    (
        "populationBalance", populationBalanceProperties, phi
    )
);
integrator->addSystem(populationBalance());
