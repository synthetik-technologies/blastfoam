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

        surfaceScalarField alphaAmaxSf(amaxSf*max(fvc::interpolate(phases[phasei]), scalar(0)));

        scalarField sumAlphaAmaxSf
        (
            fvc::surfaceSum(alphaAmaxSf)().primitiveField()
        );

        sumPhi += sumAlphaAmaxSf;

        phaseCoNums[phasei] =
            0.5*gMax(sumAlphaAmaxSf/V)*runTime.deltaTValue();
        meanPhaseCoNums[phasei] =
            0.5*(gSum(sumAlphaAmaxSf)/gSum(V))*runTime.deltaTValue();
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
    CoNum = maxCoNum;

    bool hasMassTransfer = false;
    forAll(phases, phasei)
    {
        const phaseModel& phase = phases[phasei];

        if (fluid.hasMassTransfer(phase))
        {
            hasMassTransfer = true;
            break;
        }
    }

    mDotCoNum = 0.0;
    if (hasMassTransfer)
    {
        Info<< "Maximum phase Courant numbers based on mass transfer:" << endl
            << incrIndent;

        const dimensionedScalar zeroMDot(dimDensity/dimTime, 0.0);
        forAll(phases, phasei)
        {
            const phaseModel& phase = phases[phasei];

            if (fluid.hasMassTransfer(phase))
            {
                volScalarField totalMDot
                (
                    volScalarField::New
                    (
                        IOobject::groupName("mDot", phase.name()),
                        mesh,
                        dimensionedScalar(dimDensity/dimTime, 0.0)
                    )
                );
                forAll(phases, phasej)
                {
                    const phaseModel& otherPhase = phases[phasej];
                    if (&otherPhase != &phase)
                    {
                        totalMDot += fluid.mDot(phase, otherPhase);
                    }
                }
                scalar mDotCo =
                    mag
                    (
                        gMaxMagSqr
                        (
                            (
                                totalMDot*runTime.deltaTValue()
                               /max(phase.alphaRho(), phase.residualAlphaRho())
                            )()
                        )
                    );
                Info<< indent << phase.name() << ": " << mDotCo << nl;
                mDotCoNum = max(mDotCoNum, mDotCo);
            }
        }
        Info<< endl << decrIndent;
    }
}
