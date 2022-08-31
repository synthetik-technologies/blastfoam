
// Initialise solid field pointer lists
PtrList<solidBlastThermo> thermoSolid(solidRegions.size());
PtrList<Foam::fvModels> fvModelsSolid(solidRegions.size());
PtrList<fvConstraints> fvConstraintsSolid(solidRegions.size());

// Populate solid field pointer lists
forAll(solidRegions, i)
{
    Info<< "*** Reading solid mesh thermophysical properties for region "
        << solidRegions[i].name() << nl << endl;

    Info<< "    Adding to thermoSolid\n" << endl;
    thermoSolid.set
    (
        i,
        solidBlastThermo::New
        (
            solidRegions[i],
            IOdictionary
            (
                IOobject
                (
                    "thermophysicalProperties",
                    runTime.constant(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            )
        )
    );

    Info<< "    Adding to fvModelsSolid\n" << endl;
    fvModelsSolid.set
    (
        i,
        new Foam::fvModels(solidRegions[i])
    );

    Info<< "    Adding fvConstraintsSolid\n" << endl;
    fvConstraintsSolid.set
    (
        i,
        new fvConstraints(solidRegions[i])
    );
}
