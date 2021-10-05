    const wordList solidsNames(rp["solid"]);

    PtrList<dynamicBlastFvMesh> solidRegions(solidsNames.size());

    forAll(solidsNames, i)
    {
        Info<< "Create solid mesh for region " << solidsNames[i]
            << " for time = " << runTime.timeName() << nl << endl;

        solidRegions.set
        (
            i,
            dynamicBlastFvMesh::New
            (
                IOobject
                (
                    solidsNames[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            ).ptr()
        );

        // Force calculation of geometric properties to prevent it being done
        // later in e.g. some boundary evaluation
        //(void)solidRegions[i].weights();
        //(void)solidRegions[i].deltaCoeffs();
    }
