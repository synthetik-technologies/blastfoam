//- Update internal error fields of all regions
forAll(fluidRegions, i)
{
    fluidRegions[i].updateError();
}
forAll(solidRegions, i)
{
    solidRegions[i].updateError();
}

//- Update error field boundaries
//  Done after all region errors have been updated to make sure
//  all fields are up to date
forAll(fluidRegions, i)
{
    fluidRegions[i].updateErrorBoundaries();
}
forAll(solidRegions, i)
{
    solidRegions[i].updateErrorBoundaries();
}

//- Update Meshes and check if balancing has occurred
forAll(fluidRegions, i)
{
    Info<<"Refining " << fluidRegions[i].name() << endl;
    fluidRegions[i].refine(false);
}
forAll(solidRegions, i)
{
    Info<<"Refining " << solidRegions[i].name() << endl;
    solidRegions[i].refine(false);
}

