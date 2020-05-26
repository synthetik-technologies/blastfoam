// Initialise solid field pointer lists
PtrList<volScalarField> solidRhos(solidRegions.size());
PtrList<volScalarField> solidPs(solidRegions.size());
PtrList<volScalarField> solidTs(solidRegions.size());
PtrList<volScalarField> solides(solidRegions.size());
PtrList<coordinateSystem> coordinates(solidRegions.size());
PtrList<solidThermoModel> thermos(solidRegions.size());
PtrList<radiationModel> radiations(solidRegions.size());
PtrList<fv::options> solidHeatSources(solidRegions.size());
PtrList<volScalarField> betavSolid(solidRegions.size());
PtrList<volSymmTensorField> aniAlphas(solidRegions.size());
PtrList<volScalarField> solidKappas(solidRegions.size());

// Populate solid field pointer lists
forAll(solidRegions, i)
{
    Info<< "*** Reading solid mesh thermophysical properties for region "
        << solidRegions[i].name() << nl << endl;

    solidRhos.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                solidRegions[i]
            ),
            solidRegions[i],
            dimensionedScalar("rho", dimDensity, 0.0)
        )
    );
    solidPs.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );
    solidTs.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );
    solides.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "e",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i],
            dimensionedScalar("e", sqr(dimVelocity), 0.0),
            solidThermoModel::eBoundaryTypes(solidTs[i]),
            solidThermoModel::eBoundaryBaseTypes(solidTs[i])
        )
    );
    solidKappas.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "kappa",
                runTime.timeName(),
                solidRegions[i]
            ),
            solidRegions[i],
            dimensionedScalar("zero", dimensionSet(1, 1, -3, -1, 0, 0, 0), small)
        )
    );

    Info<< "Reading thermophysicalProperties dictionary\n" << endl;
    IOdictionary thermophysicalProperties
    (
        IOobject
        (
            "thermophysicalProperties",
            runTime.constant(),
            solidRegions[i],
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info<< "    Adding to thermos\n" << endl;
    thermos.set
    (
        i,
        solidThermoModel::New
        (
            word::null,
            solidPs[i],
            solidRhos[i],
            solides[i],
            solidTs[i],
            thermophysicalProperties,
            true
        )
    );

    Info<< "    Adding to radiations\n" << endl;
    radiations.set(i, radiationModel::New(thermos[i].T()));

    Info<< "    Adding fvOptions\n" << endl;
    solidHeatSources.set
    (
        i,
        new fv::options(solidRegions[i])
    );

    if (!thermos[i].isotropic())
    {
        Info<< "    Adding coordinateSystems\n" << endl;
        coordinates.set
        (
            i,
            coordinateSystem::New(solidRegions[i], thermophysicalProperties)
        );

        tmp<volVectorField> tkappaByCp =
            thermos[i].Kappa()/thermos[i].Cp();

        aniAlphas.set
        (
            i,
            new volSymmTensorField
            (
                IOobject
                (
                    "Anialpha",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                solidRegions[i],
                dimensionedSymmTensor
                (
                    "zero",
                    tkappaByCp().dimensions(),
                    Zero
                ),
                zeroGradientFvPatchSymmTensorField::typeName
            )
        );

        aniAlphas[i].primitiveFieldRef() =
            coordinates[i].R().transformVector(tkappaByCp());
        aniAlphas[i].correctBoundaryConditions();

    }
    else
    {
        solidKappas[i] =
            solidRegions[i].lookupObject<solidThermoModel>("solidThermo").kappa();
    }

    IOobject betavSolidIO
    (
        "betavSolid",
        runTime.timeName(),
        solidRegions[i],
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (betavSolidIO.typeHeaderOk<volScalarField>(true))
    {
        betavSolid.set
        (
            i,
            new volScalarField(betavSolidIO, solidRegions[i])
        );
    }
    else
    {
        betavSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "betavSolid",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                solidRegions[i],
                dimensionedScalar(dimless, scalar(1))
            )
        );
    }
}
