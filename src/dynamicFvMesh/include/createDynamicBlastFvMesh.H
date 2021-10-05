    Info<< "Create mesh for time = "
        << runTime.timeName() << nl << endl;

    autoPtr<dynamicBlastFvMesh> meshPtr
    (
        dynamicBlastFvMesh::New
        (
            IOobject
            (
                dynamicBlastFvMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );

    dynamicBlastFvMesh& mesh = meshPtr();
