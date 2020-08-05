solidRegions[i].update();

if (thermalStress)
{
    volScalarField& T = thermo.T();
    fvScalarMatrix TEqn
    (
        fvm::ddt(T) == fvm::laplacian(DT, T) + fvOptions(T)
    );

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);
}

{
    fvVectorMatrix DEqn
    (
        fvm::d2dt2(D)
     ==
        fvm::laplacian(2*mu + lambda, D, "laplacian(DD,D)")
      + divSigmaExp
      + fvOptions.d2dt2(D)
    );

    if (thermalStress)
    {
        const volScalarField& T = thermo.T();
        DEqn += fvc::grad(threeKalpha*T);
    }

    fvOptions.constrain(DEqn);

    if (!compactNormalStress)
    {
        divSigmaExp = fvc::div(DEqn.flux());
    }

    DEqn.solve();
}

{
    volTensorField gradD(fvc::grad(D));
    sigmaD = mu*twoSymm(gradD) + (lambda*I)*tr(gradD);

    if (compactNormalStress)
    {
        divSigmaExp = fvc::div
        (
            sigmaD - (2*mu + lambda)*gradD,
            "div(sigmaD)"
        );
    }
    else
    {
        divSigmaExp += fvc::div(sigmaD);
    }
}

thermo.correct();

Info<< "max(T): " << max(thermo.T()).value()
    << ", min(T): " << min(thermo.T()).value() << endl;

if (runTime.writeTime())
{
    volSymmTensorField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()*sigmaD
    );

    if (thermalStress)
    {
        sigma = sigma - I*(thermo.rho()*threeKalpha*thermo.T());
    }

    volScalarField sigmaEq
    (
        IOobject
        (
            "sigmaEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((3.0/2.0)*magSqr(dev(sigma)))
    );

    Info<< "Max sigmaEq = " << max(sigmaEq).value()
        << endl;

    runTime.write();
}
