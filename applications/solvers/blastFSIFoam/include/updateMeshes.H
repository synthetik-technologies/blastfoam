//- Update Meshes and check if mesh moved
forAll(fluidRegions, i)
{
    Info<<"Updating mesh for " << fluidRegions[i].name() << endl;
    fluidRegions[i].update();
}
forAll(solidRegions, i)
{
    Info<<"Updating mesh for " << solidRegions[i].name() << endl;
    solidRegions[i].update();
}

