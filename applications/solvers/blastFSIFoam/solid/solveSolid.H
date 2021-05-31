
int iCorr = 0;
scalar initialResidual = 0;

tmp<volScalarField> threeKalphaTmp;

do
{
    if (thermalStress)
    {
        if (planeStress)
        {
            threeKalphaTmp = E/(1.0 - nu)*alphas[i];
        }
        else
        {
            threeKalphaTmp = E/(1.0 - 2.0*nu)*alphas[i];
        }
        volScalarField DT("DT", thermo.kappa()/thermo.Cp());

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
            DEqn += fvc::grad(threeKalphaTmp()*T);
        }

        fvOptions.constrain(DEqn);

        initialResidual = DEqn.solve().max().initialResidual();

        if (!compactNormalStress)
        {
            divSigmaExp = fvc::div(DEqn.flux());
        }
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

    Info<< "Solid Iter " << iCorr << ": residual=" << initialResidual << endl;

} while (initialResidual > convergenceTolerance && ++iCorr < nCorr);

thermo.correct();

Info<< "max(T): " << max(thermo.T()).value()
    << ", min(T): " << min(thermo.T()).value() << endl;

sigmas[i] = thermo.rho()*sigmaD;

if (thermalStress)
{
    sigmas[i] = sigmas[i] - I*(thermo.rho()*threeKalphaTmp*thermo.T());
}

sigmaEqs[i] = sqrt((3.0/2.0)*magSqr(dev(sigmas[i])));

volPointInterpolation::New(mesh).interpolate
(
    D,
    pointD
);

Info<< "Max sigmaEq = " << max(sigmaEqs[i]).value()
    << endl;

