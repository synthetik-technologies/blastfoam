#include "dictionary.H"
#include "argList.H"
#include "IntegratorsFwd.H"
#include "createEquations.H"

using namespace Foam;


createUnivariateEquation2
(
    testEqn1,
    scalar,
    0.0, 2.5,
    Foam::cos(x) - Foam::pow3(x),
    -Foam::sin(x) - 3.0*Foam::sqr(x),
    -Foam::cos(x) - 6.0*x
);

createUnivariateEquation2
(
    testEqn2,
    vector,
    0.0, 2.5,
    cmptPow(vector(x, x, x), vector::one*3.0),
    3.0*cmptPow(vector(x, x, x), vector::one*2.0),
    6.0*vector(x, x, x)
);

int main(int argc, char *argv[])
{

    testEqn1 eqn1;
    testEqn2 eqn2;
    dictionary dict;

    labelList nIntervals({1, 5, 10, 25, 50});

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
        autoPtr<scalarIntegrator> integrator
        (
            scalarIntegrator::New(eqn1, dict)
        );

        forAll(nIntervals, inti)
        {
            integrator->setNIntervals(nIntervals[inti]);
            scalar est = integrator->integrate(0, 2.5, 0);
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
    const vector ans2(9.76563*vector::one);
    forAll(intMethods2, i)
    {
        dict.set("integrator", intMethods2[i]);
        autoPtr<vectorIntegrator> integrator
        (
            vectorIntegrator::New(eqn2, dict)
        );

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
