{
    Info<< "Maximum phase Courant numbers based on eigenvalues:" << endl;

    scalarField sumPhi(mesh.V().field().size(), 0.0);
    forAll(phases, phasei)
    {
        scalarField phii
        (
            fvc::surfaceSum
            (
                mag(phases[phasei].phi())
              + fvc::interpolate(phases[phasei].speedOfSound())*mesh.magSf()
            )().primitiveField()
        );
        sumPhi = max(sumPhi, phii);
        scalar CoNumi(0.5*gMax(phii/mesh.V().field())*runTime.deltaTValue());

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
