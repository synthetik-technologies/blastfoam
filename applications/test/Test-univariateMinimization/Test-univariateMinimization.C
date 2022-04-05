#include "dictionary.H"
#include "univariateMinimizationScheme.H"
#include "argList.H"
#include "IntegratorsFwd.H"
#include "createEquations.H"

using namespace Foam;


createEquation2
(
    testEqn1,
    1.0, 3.0,
    mag(x - 2.0) + sqr(x - 1.0),
    (x - 2.0)/max(mag(x - 2.0), small) + 2.0*(x - 1.0),
    2.0
);

createEquation2
(
    testEqn2,
    1.0, 5.0,
    sqr(x - 2.0),
    2.0*(x - 2.0),
    2.0
);



int main(int argc, char *argv[])
{
    minimizationScheme::debug = 0;
    univariateMinimizationScheme::debug = 0;

    PtrList<equation> uniEqns(2);
    uniEqns.set(0, new testEqn1());
    uniEqns.set(1, new testEqn2());

    dictionary dict;

    Info<< nl << "Univariate minimization" << endl;
    wordList methods
    (
        minimizationScheme::dictionaryUnivariateConstructorTablePtr_->toc()
    );
    forAll(uniEqns, eqni)
    {
        Info<< "Solving " << uniEqns[eqni].printfx() << endl;

        const equation& eqn = uniEqns[eqni];
        forAll(methods, i)
        {
            dict.set("solver", methods[i]);
            autoPtr<minimizationScheme> solver
            (
                minimizationScheme::New(eqn, dict)
            );
            Info<< "    " << methods[i] <<": "
                << " min=" << solver->solve()
                << ", nSteps=" << solver->nSteps()
                << ", error=" << solver->errors() << endl;
        }
        Info<<nl;
    }

    Info<< "done" << endl;

    return 0;
}
