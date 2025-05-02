{
    const scalarField& V = mesh.V();

    scalarField sumPhi(V.size(), 0.0);
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

        sumPhi +=
            fvc::surfaceSum
            (
                amaxSf*fvc::interpolate(phases[phasei])
            )().primitiveField();

        phaseCoNums[phasei] =
            0.5*gMax(sumAmaxSf/V)*runTime.deltaTValue();
        meanPhaseCoNums[phasei] =
            0.5*(gSum(amaxSf)/gSum(V))*runTime.deltaTValue();
    }

    scalar maxCoNum = 0.5*gMax(sumPhi/V)*runTime.deltaTValue();
    scalar meanCoNum
    (
        0.5*(gSum(sumPhi)/gSum(V))*runTime.deltaTValue()
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
}
