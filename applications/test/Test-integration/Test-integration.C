#include "dictionary.H"
#include "argList.H"
#include "IntegratorsFwd.H"
#include "MultivariateIntegratorsFwd.H"
#include "createEquations.H"

using namespace Foam;


createNamedEquation0
(
    testEqn1,
    "f(x) = cos(x) - x^3",
    scalar,
    0.0, 2.5,
    Foam::cos(x) - Foam::pow3(x)
);

createNamedEquation0
(
    testEqn2,
    "f(X) = <x^(1.257), 3.2*x*sin(x), exp(-x)>",
    vector,
    0.0, 2.5,
    vector(Foam::pow(x, 0.1257), 3.2*x*Foam::sin(x), Foam::exp(-x))
);


createNamedUnivariateEquation0
(
    testEqn3,
    "f(x, y, z) = y*z(x + exp(y))^2 - cos(z)",
    scalar,
    scalarList(3, -great), scalarList(3, great),
    sqr(x[0] + Foam::exp(x[1]))*x[1]*x[2] - Foam::cos(x[2]);
);


int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    testEqn1 eqn1;
    testEqn2 eqn2;
    dictionary dict;

    labelList nIntervals({1, 5, 10, 25, 100, 500});

    wordList intMethods1
    (
        scalarIntegrator::dictionaryConstructorTablePtr_->toc()
    );

    Info<< "Integration of " << eqn1.eqnString() << " in (0, 2.5)" << endl;
    const scalar ans1 = -9.167152855896045;
    forAll(intMethods1, i)
    {
        dict.set("integrator", intMethods1[i]);
        dict.set
        (
            "maxSplits",
            round(Foam::log(scalar(max(nIntervals)))/Foam::log(2.0))
        );
        dict.set("nNodes", 8);
        autoPtr<scalarIntegrator> integrator
        (
            scalarIntegrator::New(eqn1, dict)
        );
        scalar est = integrator->integrate(0.0, 2.5, 0);
        Info<< incrIndent
            << indent << integrator->nIntervals() << " adaptive intervals:" << nl
            << incrIndent
            << indent << "estimate=" << est << nl
            << indent << "error=" << mag(est - ans1) << nl
            << indent << "nEvals=" << integrator->nEvals() << nl
            << decrIndent <<endl;
        forAll(nIntervals, inti)
        {
            integrator->setNIntervals(nIntervals[inti]);
            est = integrator->integrate(0.0, 2.5, 0);
            Info<< indent
                << integrator->nIntervals() << " steps: "
                << "estimate=" << est << ", error=" << mag(est - ans1) << ", "
                << integrator->nEvals() << " evals" << endl;
        }
        decrIndent(Info);
        Info<< endl;
    }

    wordList intMethods2
    (
        vectorIntegrator::dictionaryConstructorTablePtr_->toc()
    );

    Info<< nl<< "Integration of " << eqn2.eqnString() << " in (0, 2.5)" << endl;
    const vector ans2
    (
        2.491944667998765,
        8.324259785508131,
        0.917915001376102
    );
    forAll(intMethods2, i)
    {
        dict.set("integrator", intMethods2[i]);
        dict.set
        (
            "maxSplits",
            round(Foam::log(scalar(max(nIntervals)))/Foam::log(2.0))
        );
        dict.set("nNodes", 10);
        autoPtr<vectorIntegrator> integrator
        (
            vectorIntegrator::New(eqn2, dict)
        );

        vector est = integrator->integrate(0.0, 2.5, 0);
        Info<< incrIndent
            << indent << integrator->nIntervals() << " adaptive intervals:" << nl
            << incrIndent
            << indent << "estimate=" << est << nl
            << indent << "error=" << mag(est - ans2) << nl
            << indent << "nEvals=" << integrator->nEvals() << nl
            << decrIndent << endl;
        forAll(nIntervals, inti)
        {
            integrator->setNIntervals(nIntervals[inti]);
            est = integrator->integrate(0.0, 2.5, 0);
            Info<< indent
                << integrator->nIntervals() << " steps: "
                << "estimate=" << est << ", error=" << mag(est - ans2) << ", "
                << integrator->nEvals() << " evals" << endl;
        }
        decrIndent(Info);
        Info<< endl;
    }


    Info<< nl << "Multivariate integration" << endl;
    dict.clear();
    testEqn3 eqn3;
    wordList multiIntMethods
    (
        scalarMultivariateIntegrator::dictionaryConstructorTablePtr_->toc()
    );

    Info<< "Integration of " << eqn3.eqnString() << " from (0, 0, 0) : (1, 1, 1)"
        << endl;

    scalarField x0(3, 0.0);
    scalarField x1(3, 1.0);
    const scalar ans3 = 0.790494360891768;
    forAll(multiIntMethods, i)
    {
        dict.set("integrator", multiIntMethods[i]);
        dict.set("tolerance", scalarField(3, 1e-6));
        dict.set("nNodes", labelList(3, 3));
        autoPtr<scalarMultivariateIntegrator> integrator
        (
            scalarMultivariateIntegrator::New(eqn3, dict)
        );
        scalar est = integrator->integrate(x0, x1, 0);
        Info<< incrIndent
            << indent << "estimate=" << est << nl
            << indent << "error=" << mag(est - ans3) << nl
            << indent << "nEvals=" << integrator->nEvals() << nl
            << decrIndent <<endl;
    }

    Info<< nl << "done" << endl;

    return 0;
}
