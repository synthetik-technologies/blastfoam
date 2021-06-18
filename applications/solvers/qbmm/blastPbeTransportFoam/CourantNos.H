// Compute Courant number based on node velocities
scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

volScalarField speedOfSound("speedOfSound", fluid->speedOfSound());

{
    surfaceScalarField amaxSf
    (
        mag(fvc::flux(U)) + fvc::interpolate(speedOfSound)*mesh.magSf()
    );

    scalarField sumAmaxSf
    (
        fvc::surfaceSum(amaxSf)().primitiveField()
    );


    CoNum = 0.5*gMax(sumAmaxSf/mesh.V().field())*runTime.deltaTValue();
    CoNum = max(CoNum, populationBalance->CoNum());

    meanCoNum =
        0.5*(gSum(sumAmaxSf)/gSum(mesh.V().field()))*runTime.deltaTValue();
}

Info<< "Mean and max Courant Numbers = "
    << meanCoNum << " " << CoNum << endl;

// Set max Courant number based on scheme
scalar pbeCoNum = populationBalance->realizableCo();
maxCo = min(maxCo, pbeCoNum);

Info<< "Pbe Realizable Courant Number = " << pbeCoNum << endl;

Info<< "Max Courant Number = " << CoNum << endl;
