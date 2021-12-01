//- Update Meshes and check if balancing has occurred
forAll(fluidRegions, i)
{
    timeIntegrators[i].preUpdateMesh();
    refineMesh(fluidRegions[i]);
}
forAll(solidRegions, i)
{
    refineMesh(solidRegions[i]);
}
