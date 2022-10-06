PtrList<regionSolver> regions(regionTypes.size());
forAll(regionMeshes, i)
{
    regions.set
    (
        i,
        regionSolver::New
        (
            regionTypes[regionMeshes[i].name()],
            regionMeshes[i]
        ).ptr()
    );
}