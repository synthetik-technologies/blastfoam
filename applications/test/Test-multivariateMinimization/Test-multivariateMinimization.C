#include "dictionary.H"
#include "multivariateMinimizationScheme.H"
#include "argList.H"

using namespace Foam;


class multivariateTestEqn1
:
    public scalarMultivariateEquation
{
public:
    // Constructors
    multivariateTestEqn1()
    :
        scalarMultivariateEquation
        (
            scalarList(3, -30.0),
            scalarList(3, 30.0)
        )
    {}

    //- Destructor
    virtual ~multivariateTestEqn1()
    {}

    virtual label nEqns() const
    {
        return 1;
    }
    virtual label nDerivatives() const
    {
        return 0;
    }
    virtual void f
    (
        const scalarList& x,
        const label li,
        scalarList& fx
    ) const
    {
        fx[0] = sqr((x[0] + 2.0)*Foam::sin(x[0])) + sqr(x[1] + 2.0) + sqr(x[2] + 3.0);
        // fx[0] = x[0]*sqr(x[1]);
    }
};


int main(int argc, char *argv[])
{
    PtrList<scalarMultivariateEquation> multEqns(1);
    wordList multiEqnStrs(1);
    multEqns.set(0, new multivariateTestEqn1());
    multiEqnStrs[0] =
        word("    f(x1, x2) = ((x1 + 1)*sin(x))^2 + (x2 + 2)^2 + (x3 + 3)^3\n");

    dictionary dict;
    // dict.add("maxSteps", 10);

    Info<< "Mulitvariate minimization" << endl;
    multivariateMinimizationScheme::debug = 3;
    wordList multivariateMethods
    (
        multivariateMinimizationScheme::dictionaryOneConstructorTablePtr_->toc()
    );
    forAll(multEqns, eqni)
    {
        Info<< "Solving equations: " << nl << multiEqnStrs[eqni] << endl;
        scalarList x0({2, 2, 2});
        const scalarMultivariateEquation& eqns = multEqns[eqni];
        forAll(multivariateMethods, i)
        {
            dict.set("solver", multivariateMethods[i]);
            autoPtr<multivariateMinimizationScheme> solver
            (
                multivariateMinimizationScheme::New(eqns, dict)
            );
            scalarField localMin(solver->solve(x0));
            Info<< "    Local minimums=" << localMin << endl;
        }
    }

    Info<< nl << "done" << endl;

    return 0;
}
