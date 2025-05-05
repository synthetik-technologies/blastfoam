fvTimeIntegrator integrator(mesh);

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
psiuCompressibleSystem fluid(physicalProperties, mesh);
integrator.addSystem(fluid);

psiuMulticomponentThermo& thermo = fluid.rThermo();

const ::Foam::fvModels& fvModels = integrator.models();
