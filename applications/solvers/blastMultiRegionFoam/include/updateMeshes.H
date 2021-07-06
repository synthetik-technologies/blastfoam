forAll(fluidRegions, i)
{
    fluidRegions[i].update();
}
forAll(solidRegions, i)
{
    solidRegions[i].update();
}
