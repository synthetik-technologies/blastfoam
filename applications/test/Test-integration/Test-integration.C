#include "dictionary.H"
#include "argList.H"
#include "IntegratorsFwd.H"
#include "createEquations.H"

using namespace Foam;


createUnivariateEquation0
(
    testEqn1,
    scalar,
    0.0, 2.5,
    Foam::cos(x) - Foam::pow3(x)
);

createUnivariateEquation0
(
    testEqn2,
    vector,
    0.0, 2.5,
    vector(Foam::pow(x, 0.1257), 3.2*x*Foam::sin(x), Foam::exp(-x))
);

int main(int argc, char *argv[])
{

    testEqn1 eqn1;
    testEqn2 eqn2;
    dictionary dict;

    labelList nIntervals({1, 5, 10, 25, 100});

    wordList intMethods1
    (
        scalarIntegrator::dictionaryConstructorTablePtr_->toc()
    );

    Info<< "Integration of " << eqn1.printfx() << " in (0, 2.5)"
        << endl;
    const scalar ans1 = -9.167152855896045;
    forAll(intMethods1, i)
    {
        dict.set("integrator", intMethods1[i]);
        dict.set("tolerance", 1e-9);
        dict.set("maxSteps", max(nIntervals));
        dict.set("nNodes", 2);
        autoPtr<scalarIntegrator> integrator
        (
            scalarIntegrator::New(eqn1, dict)
        );
        scalar est = integrator->integrate(0.0, 2.5, 0);
        Info<< "\t" << integrator->nIntervals() << " adaptive intervals: "
            << "estimate=" << est << ", error=" << mag(est - ans1) <<endl;
        forAll(nIntervals, inti)
        {
            integrator->setNIntervals(nIntervals[inti]);
            est = integrator->integrate(0.0, 2.5, 0);
            Info<< "\t" << integrator->nIntervals() << " steps: "
                << "estimate=" << est << ", error=" << mag(est - ans1) <<endl;
        }
        Info<< endl;
    }

    wordList intMethods2
    (
        vectorIntegrator::dictionaryConstructorTablePtr_->toc()
    );

    Info<< "Integration of " << eqn2.printfx() << " in (0, 2.5)"
        << endl;
    const vector ans2
    (
        2.491944667998765,
        8.324259785508131,
        0.917915001376102
    );
    forAll(intMethods2, i)
    {
        dict.set("integrator", intMethods2[i]);
        dict.set("tolerance", 1e-9);
        dict.set("maxSteps", max(nIntervals)*2);
        dict.set("nNodes", 5);
        autoPtr<vectorIntegrator> integrator
        (
            vectorIntegrator::New(eqn2, dict)
        );

        vector est = integrator->integrate(0.0, 2.5, 0);
        Info<< "\t" << integrator->nIntervals() << " adaptive steps: "
            << "estimate=" << est << ", error=" << mag(est - ans2) <<endl;
        forAll(nIntervals, inti)
        {
            integrator->setNIntervals(nIntervals[inti]);
            vector est = integrator->integrate(0, 2.5, 0);
            Info<< "\t" << integrator->nIntervals() << " steps: "
                << "estimate=" << est << ", error=" << mag(est - ans2) <<endl;
        }
        Info<< endl;
    }

    Info<< nl << "done" << endl;

    return 0;
}
