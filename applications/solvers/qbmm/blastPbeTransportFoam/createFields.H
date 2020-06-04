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

Info<< "Reading MieGruneisenEOSProperties dictionary\n" << endl;
IOdictionary MieGruneisenEOSProperties
(
     IOobject
     (
          "MieGruneisenEOSProperties",
          runTime.constant(),
          mesh,
          IOobject::MUST_READ,
          IOobject::NO_WRITE
     )
);

autoPtr<phaseCompressibleSystem> fluid
(
    phaseCompressibleSystem::New(mesh, MieGruneisenEOSProperties)
);
autoPtr<timeIntegrator> integrator(timeIntegrator::New(mesh));
integrator->addSystem(fluid());

volScalarField& p = fluid->p();
const volVectorField& U = fluid->U();
const volScalarField& rho = fluid->rho();
const surfaceScalarField& phi = fluid->phi();
fluid->update();

autoPtr<errorEstimator> error;
if (!isA<staticFvMesh>(mesh))
{
    error =
        errorEstimator::New
        (
            mesh,
            IOdictionary
            (
                IOobject
                (
                    "dynamicMeshDict",
                    runTime.constant(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );
}



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


autoPtr<populationBalanceModel> populationBalance
(
    populationBalanceModel::New
    (
        "populationBalance", populationBalanceProperties, phi
    )
);
