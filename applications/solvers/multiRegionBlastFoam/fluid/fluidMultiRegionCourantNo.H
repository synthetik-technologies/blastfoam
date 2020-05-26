scalar CoNum = -great;

forAll(fluidRegions, regionI)
{
    volScalarField speedOfSound("speedOfSound", fluids[regionI].speedOfSound());

    surfaceScalarField amaxSf
    (
        mag(fluids[regionI].phi())
      + fvc::interpolate(speedOfSound)*fluidRegions[regionI].magSf()
    );
    scalarField sumAmaxSf
    (
        fvc::surfaceSum(amaxSf)().primitiveField()
    );

    CoNum =
        max
        (
            CoNum,
            0.5
           *gMax(sumAmaxSf/fluidRegions[regionI].V().field())
           *runTime.deltaTValue()
        );
}
Info<< "Max Courant Numbers = " << CoNum << endl;
