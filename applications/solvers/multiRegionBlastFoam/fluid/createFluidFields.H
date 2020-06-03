// Initialise fluid field pointer lists
PtrList<phaseCompressibleSystem> fluids(fluidRegions.size());
PtrList<timeIntegrator> integrators(fluidRegions.size());
PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
PtrList<errorEstimator> errors(fluidRegions.size());
PtrList<volScalarField> kappas(fluidRegions.size());

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    Info<< "    Adding to gFluid\n" << endl;
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

    kappas.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "kappa",
                runTime.timeName(),
                fluidRegions[i]
            ),
            fluidRegions[i],
            dimensionedScalar("zero", dimensionSet(1, 1, -3, -1, 0, 0, 0), small)
        )
    );

    Info<< "    Adding to fluids\n" << endl;
    fluids.set
    (
        i,
        phaseCompressibleSystem::New(fluidRegions[i], phaseProperties).ptr()
    );
    kappas[i] = fluidRegions[i].lookupObject<fluidThermoModel>("fluidThermo").kappa();
    integrators.set
    (
        i,
        timeIntegrator::New(fluidRegions[i]).ptr()
    );
    integrators[i].addSystem(fluids[i]);



    fluids[i].update();



//     Info<< "    Adding to errorEstimators\n" << endl;
//     errors.set
//     (
//         i,
//         errorEstimator::New
//         (
//             fluidRegions[i],
//             IOdictionary
//             (
//                 IOobject
//                 (
//                     "dynamicMeshDict",
//                     runTime.constant(),
//                     fluidRegions[i],
//                     IOobject::MUST_READ,
//                     IOobject::NO_WRITE
//                 )
//             )
//         ).ptr()
//     );
}
