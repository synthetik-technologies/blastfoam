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

HashTable<word> regionTypes;
forAllConstIter(HashTable<wordList>, regionProperties, iter)
{
    const wordList& regionNames(iter());
    forAll(regionNames, i)
    {
        regionTypes.insert(regionNames[i], iter.key());
    }
}

PtrList<dynamicFvMesh> regionMeshes(regionTypes.size());
{
    label regioni = 0;
    forAllConstIter(HashTable<word>, regionTypes, iter)
    {
        Info<<"Creating mesh for region " << iter.key() << endl;
        regionMeshes.set
        (
            regioni++,
            dynamicFvMesh::New
            (
                IOobject
                (
                    iter.key(),
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            ).ptr()
        );
    }
}