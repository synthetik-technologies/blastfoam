#include "dictionary.H"
#include "univariateRootSolver.H"
#include "argList.H"
#include "createEquations.H"

using namespace Foam;


createEquation2
(
    testEqn1,
    scalar,
    0.01, 1.0,
    Foam::cos(x) - Foam::pow3(x),
    -Foam::sin(x) - 3.0*Foam::sqr(x),
    -Foam::cos(x) - 6.0*x
);

createEquation2
(
    testEqn2,
    scalar,
    0.01, 1.0,
    Foam::exp(x) - 10.0*x,
    Foam::exp(x) - 10.0,
    Foam::exp(x)
);


int main(int argc, char *argv[])
{

    PtrList<scalarEquation> uniEqns(2);
    uniEqns.set(0, new testEqn1());
    uniEqns.set(1, new testEqn2());

    dictionary dict;

    Info<< endl;
    Info<< "Univariate root finding" << endl;
    wordList methods
    (
        univariateRootSolver::dictionaryTwoConstructorTablePtr_->toc()
    );
    forAll(uniEqns, eqni)
    {
        Info<< "Solving " << uniEqns[eqni].name() << endl;
        const scalarEquation& eqn = uniEqns[eqni];
        forAll(methods, i)
        {
            dict.set("solver", methods[i]);
            autoPtr<univariateRootSolver> solver
            (
                univariateRootSolver::New(eqn, dict)
            );
            Info<<"    root=" << solver->solve(0.5, 0)
                << ", nSteps=" << solver->nSteps()
                << ", error=" << solver->error() << nl << endl;
        }
    }

    Info<< "done" << endl;

    return 0;
}
