    const wordList solidNames
    (
        rp.found("solid") ? rp["solid"] : wordList(0)
    );

    PtrList<dynamicFvMesh> solidRegions(solidNames.size());

    forAll(solidNames, i)
    {
        Info<< "Create solid mesh for region " << solidNames[i]
            << " for time = " << runTime.timeName() << nl << endl;

        solidRegions.set
        (
            i,
            dynamicFvMesh::New
            (
                IOobject
                (
                    solidNames[i],
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
