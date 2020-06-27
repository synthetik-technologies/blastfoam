// Initialise fluid field pointer lists
PtrList<phaseCompressibleSystem> fluids(fluidRegions.size());
PtrList<timeIntegrator> timeIntegrators(fluidRegions.size());
PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    // Read gravitational acceleration
    Info<< "\nReading g" << endl;
    gFluid.set
    (
        i,
        new uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedVector(dimAcceleration, Zero)
        )
    );

    Info<< "Reading phaseProperties dictionary\n" << endl;
    IOdictionary phaseProperties
    (
        IOobject
        (
            "phaseProperties",
            runTime.constant(),
            fluidRegions[i],
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    fluids.set
    (
        i,
        phaseCompressibleSystem::New(fluidRegions[i], phaseProperties)
    );
    timeIntegrators.set(i, timeIntegrator::New(fluidRegions[i]));
    timeIntegrators[i].addSystem(fluids[i]);
    fluids[i].update();
}
