regionProperties regionProperties(runTime);
IOdictionary regionPropertiesDict
(
    IOobject
    (
        "regionProperties",
        runTime.time().constant(),
        runTime.db(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

List<Pair<word>> regionTypes;
forAllConstIter(HashTable<wordList>, regionProperties, iter)
{
    const wordList& regionNames(iter());
    forAll(regionNames, regioni)
    {
        regionTypes.append({regionNames[regioni], iter.key()});
    }
}

PtrList<dynamicFvMesh> regionMeshes(regionTypes.size());
{
    forAll(regionTypes, regioni)
    {
        Info<<"Creating mesh for region "
            << regionTypes[regioni].first() << endl;
        regionMeshes.set
        (
            regioni,
            dynamicFvMesh::New
            (
                IOobject
                (
                    regionTypes[regioni].first(),
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            ).ptr()
        );
    }
}
