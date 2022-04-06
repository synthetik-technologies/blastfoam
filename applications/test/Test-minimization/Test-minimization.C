#include "dictionary.H"
#include "minimizationScheme.H"
#include "createEquations.H"
#include "UnivariateEquationsFwd.H"
#include "argList.H"

using namespace Foam;

createNamedUnivariateEquation0
(
    testEqn1,
    "f(x, y, z) = ((x + 2)*sin(x))^2 + (y + 2)^2 + (z + 3)^2",
    scalar,
    scalarList({2.0, -4.0, -4.0}), scalarList({4.0, 4.0, 4.0}),
    Foam::sqr((x[0] + 2.0)*Foam::sin(x[0]))
  + Foam::sqr(x[1] + 2.0)
  + Foam::sqr(x[2] + 3.0)
);

createNamedEquation1
(
    testEqn2,
    "f(x) = (x - 2)^2",
    scalar,
    -3, 3,
    sqr(x - 2.0),
    2.0*(x - 2.0)
);


int main(int argc, char *argv[])
{
    PtrList<scalarUnivariateEquation> multEqns(2);
    multEqns.set(0, new testEqn1());
    multEqns.set(1, new testEqn2());

    dictionary dict;
//     dict.add("maxSteps", 10);
    dict.add("nParticles", 1000);

    Info<< "Minimization" << endl;
    minimizationScheme::debug = 0;
    forAll(multEqns, eqni)
    {
        Info<< nl << "Solving equations: " << nl << multEqns[eqni].name() << endl;
        scalarField x0Orig(scalarList({0, 0, 0}));
        const scalarUnivariateEquation& eqns = multEqns[eqni];

        wordList multivariateMethods
        (
            eqns.nVar() == 1
          ? minimizationScheme::dictionaryUnivariateConstructorTablePtr_->toc()
          : minimizationScheme::dictionaryMultivariateConstructorTablePtr_->toc()
        );

        forAll(multivariateMethods, i)
        {
            scalarField x0(x0Orig);
            x0.resize(eqns.nVar());
            dict.set("solver", multivariateMethods[i]);
            autoPtr<minimizationScheme> solver
            (
                minimizationScheme::New(eqns, dict)
            );
            scalarField localMin(solver->solve());
            Info<< "    " << multivariateMethods[i] << ": "
                << "Local minimums = " << localMin << endl;
        }
    }

    Info<< nl << "done" << endl;

    return 0;
}
