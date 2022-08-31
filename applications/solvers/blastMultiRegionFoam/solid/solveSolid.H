label maxIter = solidRegions[i].solutionDict().lookup<label>("maxIter");
scalar tolerance = solidRegions[i].solutionDict().lookup<scalar>("tolerance");

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
        fvm::ddt(rho, e)
      + thermo.divq(e)
     ==
        fvModels.source(rho, e)
    );



    if (!lastIter)
    {
        eEqn.relax();
    }

    fvConstraints.constrain(eEqn);

    converged = eEqn.solve().initialResidual() < tolerance;

    if (!lastIter)
    {
        e.storePrevIter();
        e.relax();
    }

    fvConstraints.constrain(e);

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

Info<< "Min/max T:" << min(thermo.T()).value() << ' '
    << max(thermo.T()).value() << endl;
