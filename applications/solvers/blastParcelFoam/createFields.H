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
coupledMultiphaseCompressibleSystem fluid(mesh);
integrator->addSystem(fluid);

basicThermoCloud parcels
(
    "reactingCloud1",
    fluid.rho(),
    fluid.U(),
    g,
    fluid.thermo()
);
volScalarField theta("theta.reactingCloud1", parcels.theta());
fluid.setDispersedVolumeFraction(theta);


const volScalarField& p = fluid.p();
const volScalarField& T = fluid.T();
const surfaceScalarField& phi = fluid.phi();
fluid.update();
