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
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);

fvTimeIntegrator integrator(mesh);
autoPtr<compressibleSystem> fluidPtr(coupledCompressibleSystem::New(mesh));
compressibleSystem& fluid = fluidPtr();
coupledCompressibleSystem& coupledFluid =
    dynamicCast<coupledCompressibleSystem>(fluid);
integrator.addSystem(fluid);

parcelCloudList clouds
(
    fluid.rho(),
    fluid.U(),
    g,
    fluid.thermo()
);
volScalarField alpha
(
    IOobject
    (
        IOobject::groupName("alpha", parcelCloudList::cloudNamesName),
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    clouds.alpha()
);
coupledFluid.setDispersedVolumeFraction(alpha);


const volScalarField& p = fluid.p();
const volScalarField& T = fluid.T();
const surfaceScalarField& phi = fluid.phi();
fluid.update();


const ::Foam::fvModels& fvModels = integrator.models();
const ::Foam::fvConstraints& fvConstraints = integrator.constraints();
