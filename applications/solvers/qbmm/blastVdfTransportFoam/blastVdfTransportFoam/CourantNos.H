// Compute Courant number based on node velocities
scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

volScalarField speedOfSound("speedOfSound", fluid.speedOfSound());

{
    scalar rPbeCoNum = populationBalance->realizableCo();
    // Set max Courant number based on scheme
    maxCo = min(maxCo, rPbeCoNum);
    Info<< "Pbe Realizable Courant Number = " << rPbeCoNum << endl;

    scalar pbeCoNum = populationBalance->CoNum();
    Info<< "Max pbe Courant Number = " << pbeCoNum << endl;

    scalarField sumAmaxSf
    (
        fvc::surfaceSum
        (
            mag(phi) + fvc::interpolate(speedOfSound)*mesh.magSf()
        )().primitiveField()
    );
    CoNum = 0.5*gMax(sumAmaxSf/mesh.V().field())*runTime.deltaTValue();
    Info<< "Max fluid Courant Number = " << CoNum << endl;

    CoNum = max(CoNum, pbeCoNum);
}
