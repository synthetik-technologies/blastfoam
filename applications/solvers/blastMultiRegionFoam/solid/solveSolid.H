label maxIter = solidRegions[i].solutionDict().lookupType<label>("maxIter");
scalar tolerance = solidRegions[i].solutionDict().lookupType<scalar>("tolerance");

label iter = 0;
bool lastIter = false;
bool converged = false;
do
{
    Info<<"Iteration " << iter << endl;
    if (converged)
    {
        lastIter = true;
    }
    fvScalarMatrix eEqn
    (
        fvm::ddt(betav*rho, e)
      - (
            thermo.isotropic()
          ? fvm::laplacian(betav*thermo.alpha(), e, "laplacian(alpha,e)")
          : fvm::laplacian(betav*taniAlpha(), e, "laplacian(alpha,e)")
        )
//      ==
//         fvOptions(rho, e)
    );

    if (!lastIter)
    {
        e.storePrevIter();

        e.relax();
        eEqn.relax();
    }

//     fvOptions.constrain(eEqn);

    converged = eEqn.solve().initialResidual() < tolerance;

//     fvOptions.correct(e);

} while (++iter < maxIter && !lastIter);
if (converged)
{
    Info<< "Converged in " << iter << " iterations" << endl;
}
else
{
    Info<< "Did not converge in " << iter << " iterations" << endl;
}

thermo.correct();

Info<< "max(T): " << max(T).value()
    << ", min(T): " << min(T).value() << endl;
