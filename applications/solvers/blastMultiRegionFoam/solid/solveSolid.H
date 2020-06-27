{
    fvScalarMatrix eEqn
    (
        fvm::ddt(betav*rho, e)
        - (
            thermo.isotropic()
            ? fvm::laplacian(betav*thermo.alpha(), e, "laplacian(alpha,e)")
            : fvm::laplacian(betav*taniAlpha(), e, "laplacian(alpha,e)")
        )
        ==
        fvOptions(rho, e)
    );

    eEqn.relax();

    fvOptions.constrain(eEqn);

    eEqn.solve();

    fvOptions.correct(e);
}

thermo.correct();

Info<< "max(T): " << max(thermo.T()).value()
    << ", min(T): " << min(thermo.T()).value() << endl;
