PtrList<regionSolver> regions(regionMeshes.size());
forAll(regionMeshes, regioni)
{
    regions.set
    (
        regioni,
        regionSolver::New
        (
            regionTypes[regioni].second(),
            regionMeshes[regioni]
        ).ptr()
    );
}
