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
autoPtr<timeIntegrator> integrator(timeIntegrator::New(fluid()));

volScalarField& p = fluid->p();
const volVectorField& U = fluid->U();
const volScalarField& rho = fluid->rho();
const surfaceScalarField& phi = fluid->phi();
fluid->update();

volScalarField error
(
    IOobject
    (
        "error",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
    0.0,
    wordList(p.boundaryField().size(), "zeroGradient")
);

volScalarField maxP
(
    IOobject
    (
        "maxP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p
);

volScalarField impulse
(
    IOobject
    (
        "impulse",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimPressure*dimTime, 0.0)
);
volScalarField maxImpulse
(
    IOobject
    (
        "maxImpulse",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    impulse
);

List<vector> stopPoints;
scalar stopEpsilon = 0.0;
Switch stopActive = false;
if (MieGruneisenEOSProperties.found("stop"))
{
    dictionary stopDict(MieGruneisenEOSProperties.subDict("stop"));

    /// stop simluation when pertubation sare detected ...
    stopPoints = stopDict.lookupType<List<vector>>("points");
    stopEpsilon = stopDict.lookupOrDefault("epsilon", 1e-3);
    stopActive = stopDict.lookupOrDefault("active", true);

    Info << "Stop points detected:" << endl;
    forAll(stopPoints, pointI)
    {
        Info << stopPoints[pointI] << endl;
        Info << "Corresponding cell id: " << mesh.findNearestCell(stopPoints[pointI]) << endl;
    }

    Info << "Ending simulation when mag(U) > " << stopEpsilon << " " << endl;
}



