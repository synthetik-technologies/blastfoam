{
    Info<< "Maximum phase Courant numbers based on eigenvalues:" << endl;

    scalarField sumPhi(mesh.V().field().size(), 0.0);
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

        sumPhi = max(sumPhi, sumAmaxSf);
        scalar CoNumi(0.5*gMax(sumAmaxSf/mesh.V().field())*runTime.deltaTValue());

        Info<< "\t" << phases[phasei].name() << ": " << CoNumi << nl;
    }
    CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();
    scalar meanCoNum
    (
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue()
    );
    Info<<"Courant number mean: " << meanCoNum
        << " max: " << CoNum << endl;
}
