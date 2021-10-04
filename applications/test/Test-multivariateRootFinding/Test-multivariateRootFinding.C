#include "dictionary.H"
#include "multivariateRootSolver.H"
#include "argList.H"

using namespace Foam;


class multivariateTestEqn1
:
    public scalarMultivariateEquation
{
public:
    // Constructors
    multivariateTestEqn1()
    {}

    //- Destructor
    virtual ~multivariateTestEqn1()
    {}

    virtual label nEqns() const
    {
        return 2;
    }
    virtual void f
    (
        const scalarField& x,
        const label li,
        scalarField& fx
    ) const
    {
        fx[0] = sqr(x[0]) + sqr(x[1]) - 4.0;
        fx[1] = sqr(x[0]) - x[1] + 1.0;
    }
    virtual void jacobian
    (
        const scalarField& x,
        const label li,
        scalarField& fx,
        scalarSquareMatrix& dfdx
    ) const
    {
        f(x, li, fx);

        dfdx(0, 0) = stabilise(2.0*x[0], small);
        dfdx(0, 1) = 2.0*x[1];
        dfdx(1, 0) = 2.0*x[0];
        dfdx(1, 1) = -1.0;
    }
};


class multivariateTestEqn2
:
    public scalarMultivariateEquation
{
public:
    // Constructors
    multivariateTestEqn2()
    {}

    //- Destructor
    virtual ~multivariateTestEqn2()
    {}

    virtual label nEqns() const
    {
        return 2;
    }
    virtual void f
    (
        const scalarField& x,
        const label li,
        scalarField& fx
    ) const
    {
        fx[0] = x[1] - sqr(x[0]) + x[0];
        fx[1] = sqr(x[0])/16.0 + sqr(x[1]) - 1.0;
    }
    virtual void jacobian
    (
        const scalarField& x,
        const label li,
        scalarField& fx,
        scalarSquareMatrix& dfdx
    ) const
    {
        f(x, li, fx);

        dfdx(0, 0) = -2.0*x[0] + 1.0;
        dfdx(0, 1) = 1.0;
        dfdx(1, 0) = 2.0*x[0]/16.0;
        dfdx(1, 1) = stabilise(2.0*x[1], small);
    }
};


int main(int argc, char *argv[])
{
    PtrList<scalarMultivariateEquation> multEqns(2);
    wordList multiEqnStrs(2);
    multEqns.set(0, new multivariateTestEqn1());
    multiEqnStrs[0] =
        word("    f1(x1, x2) = x1^2 + x2^2 - 4.0\n")
      + word("    f2(x1, x2) = x1^2 - x2 + 1.0\n");
    multEqns.set(1, new multivariateTestEqn2());
    multiEqnStrs[1] =
        word("    f1(x1, x2) = x2^2 + x1^2 + x1\n")
      + word("    f2(x1, x2) = (x1^2)/16 - x2^2 - 1.0\n");

    dictionary dict;

    Info<< "Mulitvariate root finding" << endl;
    wordList multivariateMethods
    (
        multivariateRootSolver::dictionaryOneConstructorTablePtr_->toc()
    );
    forAll(multEqns, eqni)
    {
        Info<< "Solving equations: " << nl << multiEqnStrs[eqni] << endl;
        scalarField x0(2, 0.1);
        const scalarMultivariateEquation& eqns = multEqns[eqni];
        forAll(multivariateMethods, i)
        {
            dict.set("solver", multivariateMethods[i]);
            autoPtr<multivariateRootSolver> solver
            (
                multivariateRootSolver::New(eqns, dict)
            );
            Info<< "    roots=" << solver->solve(x0, 0)
                << ", nSteps=" << solver->nSteps()
                << ", error=" << solver->error() << nl << endl;
        }
    }

    Info<< "done" << endl;

    return 0;
}
