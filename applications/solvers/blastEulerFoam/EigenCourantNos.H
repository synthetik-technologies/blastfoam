{
    scalarField sumPhi(mesh.V().field().size(), 0.0);
    scalarList phaseCoNums(phases.size(), 0.0);
    scalarList meanPhaseCoNums(phases.size(), 0.0);
    forAll(phases, phasei)
    {
        surfaceScalarField amaxSf
        (
            fvc::interpolate(phases[phasei].speedOfSound())*mesh.magSf()
        );
        // Remove wave speed from wedge boundaries
        forAll(amaxSf.boundaryField(), patchi)
        {
            if (isA<wedgeFvPatch>(mesh.boundary()[patchi]))
            {
                amaxSf.boundaryFieldRef() = Zero;
            }
        }
        amaxSf += mag(fvc::flux(phases[phasei].U()));

        scalarField sumAmaxSf
        (
            fvc::surfaceSum(amaxSf)().primitiveField()
        );

        sumPhi += fvc::interpolate(phases[phasei])*sumAmaxSf;

        phaseCoNums[phasei] =
            0.5*gMax(sumAmaxSf/mesh.V().field())*runTime.deltaTValue();
        meanPhaseCoNums[phasei] =
            0.5*(gSum(amaxSf)/gSum(mesh.V().field()))*runTime.deltaTValue();
    }

    scalar maxCoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();
    scalar meanCoNum
    (
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue()
    );

    Info<< "Courant number: mean = " << meanCoNum
        << ", max = " << maxCoNum << endl;

    Info<< "Phase Courant numbers based on eigenvalues:"
        << incrIndent << endl;

    forAll(phases, phasei)
    {
        Info<< indent << phases[phasei].name() << ": "
            << "mean = " << meanPhaseCoNums[phasei]
            << ", max = " << phaseCoNums[phasei] << nl;
    }
    Info<< endl << decrIndent;
    CoNum = max(phaseCoNums);//max(max(phaseCoNums), maxCoNum);
    Info<<"CoNum "<<CoNum<<endl;
}
