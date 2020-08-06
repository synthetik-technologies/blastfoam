// Initialise solid field pointer lists

// General fields
PtrList<volVectorField> Ds(solidRegions.size());
PtrList<fv::options> solidSources(solidRegions.size());
PtrList<volVectorField> divSigmaExps(solidRegions.size());
PtrList<volSymmTensorField> sigmaDs(solidRegions.size());
PtrList<volSymmTensorField> sigmas(solidRegions.size());
PtrList<volScalarField> sigmaEqs(solidRegions.size());

// Mechanical fields
PtrList<IOdictionary> mechanicalProperties(solidRegions.size());
PtrList<volScalarField> Es(solidRegions.size());
PtrList<volScalarField> nus(solidRegions.size());
PtrList<volScalarField> mus(solidRegions.size());
PtrList<volScalarField> lambdas(solidRegions.size());
List<Switch> planeStresses(solidRegions.size(), false);
List<Switch> compactNormalStresses(solidRegions.size(), false);


// Thermal fields
PtrList<IOdictionary> thermophysicalProperties(solidRegions.size());
PtrList<solidThermoModel> solidThermos(solidRegions.size());
// PtrList<radiationModel> radiations(solidRegions.size());
PtrList<volScalarField> alphas(solidRegions.size());
List<Switch> thermalStresses(solidRegions.size(), false);


// Populate solid field pointer lists
forAll(solidRegions, i)
{
    Info<< "*** Reading solid mesh thermophysical properties for region "
    << solidRegions[i].name() << nl << endl;

    thermophysicalProperties.set
    (
        i,
        new IOdictionary
        (
            IOobject
            (
                "thermophysicalProperties",
                runTime.constant(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    Info<< "    Adding to thermos\n" << endl;
    solidThermos.set
    (
        i,
        solidThermoModel::New
        (
            word::null,
            solidRegions[i],
            thermophysicalProperties[i],
            true
        )
    );

    // Info<< "    Adding to radiations\n" << endl;
    // radiations.set(i, radiationModel::New(solidThermos[i].T()));

    #include "readMechanicalProperties.H"
    #include "readThermalProperties.H"

    Ds.set
    (
        i,
        new volVectorField
        (
             IOobject
            (
                "D",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );

    sigmaDs.set
    (
        i,
        new volSymmTensorField
        (
             IOobject
            (
                "sigmaD",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mus[i]*twoSymm(fvc::grad(Ds[i])) + lambdas[i]*(I*tr(fvc::grad(Ds[i])))
        )
    );

    divSigmaExps.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "divSigmaExp",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::div(sigmaDs[i])
        )
    );

    sigmas.set
    (
        i,
        new volSymmTensorField
        (
             IOobject
            (
                "sigma",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            solidThermos[i].rho()*sigmaDs[i]
        )
    );

    sigmaEqs.set
    (
        i,
        new volScalarField
        (
             IOobject
            (
                "sigmaEq",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sqrt((3.0/2.0)*magSqr(dev(sigmas[i])))
        )
    );

    const dictionary& stressControl =
        solidRegions[i].solutionDict().subDict("stressAnalysis");

    compactNormalStresses[i] =
        stressControl.lookupType<Switch>("compactNormalStress");

    if (compactNormalStresses[i])
    {
        divSigmaExps[i] -=
            fvc::laplacian(2*mus[i] + lambdas[i], Ds[i], "laplacian(DD,D)");
    }
    else
    {
        divSigmaExps[i] -=
            fvc::div((2*mus[i] + lambdas[i])*fvc::grad(Ds[i]), "div(sigmaD)");
    }

    solidRegions[i].setFluxRequired(Ds[i].name());

    Info<< "    Adding fvOptions\n" << endl;
    solidSources.set
    (
        i,
        new fv::options(solidRegions[i])
    );
}
