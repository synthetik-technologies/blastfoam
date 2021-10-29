{
    volScalarField& he = thermo.he();

    fvScalarMatrix deltaEEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div(fvc::absolute(phi, rho, U), p/rho)
          : -dpdt
        )
      + thermophysicalTransport->divq(he)
     ==
        fvModels.source(rho, he)
    );

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he)
      + deltaEEqn
    );

    if (IBM.thermalForcingNeeded())
    {
        EEqn -=
            IBM.forcing
            (
                he,
                rho,
                (rho.oldTime()*he.oldTime())(),
                (deltaEEqn & he)(),
                runTime.deltaT()
            );
    }
    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

    thermo.correct();
}
