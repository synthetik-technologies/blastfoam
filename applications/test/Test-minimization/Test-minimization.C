#include "dictionary.H"
#include "minimizationScheme.H"
#include "createEquations.H"
#include "EquationsFwd.H"
#include "UnivariateEquationsFwd.H"
#include "argList.H"

using namespace Foam;

createNamedEquation2
(
    testEqn1,
    "f(x) = |x - 2| + (x - 1)^2",
    scalar,
    1.0, 3.0,
    mag(x - 2.0) + sqr(x - 1.0),
    (x - 2.0)/max(mag(x - 2.0), small) + 2.0*(x - 1.0),
    2.0
);

createNamedEquation2
(
    testEqn2,
    "f(x) = (x - 2)^2",
    scalar,
    1.0, 4.0,
    sqr(x - 2.0),
    2.0*(x - 2.0),
    2.0
);

createNamedUnivariateEquation0
(
    testEqn3,
    "f(x, y, z) = ((x + 2)*sin(x))^2 + (y + 2)^2 + (z + 3)^2",
    scalar,
    scalarList({2.0, -4.0, -4.0}), scalarList({4.0, 4.1, 4.2}),
    Foam::sqr((x[0] + 2.0)*Foam::sin(x[0]))
  + Foam::sqr(x[1] + 2.0)
  + Foam::sqr(x[2] + 3.0)
);


int main(int argc, char *argv[])
{
    minimizationScheme::debug = 2;

    PtrList<scalarUnivariateEquation> eqns(3);
    eqns.set(0, new testEqn1());
    eqns.set(1, new testEqn2());
    eqns.set(2, new testEqn3());

    dictionary dict;
    dict.add("cLocal", 0.3);
    dict.add("cGlobal", 0.1);
    dict.add("vWeight", 0.8);
    dict.add("nParticles", 1000);

    dict.add("maxSteps", 200);


    Info<< nl << "Univariate minimization" << endl;
    wordList methods
    (
        minimizationScheme::dictionaryUnivariateConstructorTablePtr_->toc()
    );

    forAll(eqns, eqni)
    {
        const scalarUnivariateEquation& eqn = eqns[eqni];
        Info<< nl << "Solving equations: " << nl << eqn.name() << nl
            << "Bounds: " << eqn.lowerLimits() << ", "
            << eqn.upperLimits()
            << nl << endl;

        wordList methods
        (
            eqn.nVar() == 1
          ? minimizationScheme::dictionaryUnivariateConstructorTablePtr_->toc()
          : minimizationScheme::dictionaryMultivariateConstructorTablePtr_->toc()
        );

        forAll(methods, i)
        {
            dict.set("solver", methods[i]);
            autoPtr<minimizationScheme> solver
            (
                minimizationScheme::New(eqn, dict)
            );
            scalarField localMin(solver->solve());
            Info<< endl;
        }
        Info<<nl;
    }

    Info<< nl << "done" << endl;

    return 0;
}
