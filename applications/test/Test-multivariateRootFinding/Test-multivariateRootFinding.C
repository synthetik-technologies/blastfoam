#include "dictionary.H"
#include "rootSolver.H"
#include "argList.H"

using namespace Foam;


class testEqn1
:
    public ScalarMultivariateEquation
{
public:
    // Constructors
    testEqn1()
    :
        ScalarMultivariateEquation
        (
            word("f1(x1, x2) = x1^2 + x2^2 - 4.0\n")
          + word("f2(x1, x2) = x1^2 - x2 + 1.0\n"),
            2,
            {0.0, 0.0},
            {2.0, 2.0}
        )
    {}

    //- Destructor
    virtual ~testEqn1()
    {}

    virtual label nEqns() const
    {
        return 2;
    }
    virtual label nDerivatives() const
    {
        return 1;
    }
    virtual void FX
    (
        const scalarList& x,
        const label li,
        scalarList& fx
    ) const
    {
        fx[0] = sqr(x[0]) + sqr(x[1]) - 4.0;
        fx[1] = sqr(x[0]) - x[1] + 1.0;
    }
    // virtual void jacobian
    // (
    //     const scalarList& x,
    //     const label li,
    //     scalarList& fx,
    //     scalarRectangularMatrix& dfdx
    // ) const
    // {
    //     f(x, li, fx);
    //     dfdx(0, 0) = stabilise(2.0*x[0], small);
    //     dfdx(0, 1) = 2.0*x[1];
    //     dfdx(1, 0) = 2.0*x[0];
    //     dfdx(1, 1) = -1.0;
    // }
};


class testEqn2
:
    public ScalarMultivariateEquation
{
public:
    // Constructors
    testEqn2()
    :
        ScalarMultivariateEquation
        (
            word("f1(x1, x2) = x2^2 + x1^2 + x1\n")
          + word("f2(x1, x2) = (x1^2)/16 - x2^2 - 1.0\n"),
            2,
            {-10.0, -10.0},
            {10.0, 10.0}
        )
    {}

    //- Destructor
    virtual ~testEqn2()
    {}

    virtual label nEqns() const
    {
        return 2;
    }
    virtual label nDerivatives() const
    {
        return 1;
    }
    virtual void FX
    (
        const scalarList& x,
        const label li,
        scalarList& fx
    ) const
    {
        fx[0] = x[1] - sqr(x[0]) + x[0];
        fx[1] = sqr(x[0])/16.0 + sqr(x[1]) - 1.0;
    }
    virtual void jacobian
    (
        const scalarList& x,
        const label li,
        scalarList& fx,
        RectangularMatrix<scalar>& J
    ) const
    {
        FX(x, li, fx);

        J(0, 0) = stabilise(-2.0*x[0] + 1.0, small);
        J(0, 1) = 1.0;
        J(1, 0) = 2.0*x[0]/16.0;
        J(1, 1) = stabilise(2.0*x[1], small);
    }
};


int main(int argc, char *argv[])
{
    PtrList<multivariateEquation<scalar>> multEqns(2);
    multEqns.set(0, new testEqn1());
    multEqns.set(1, new testEqn2());

    dictionary dict;

    Info<< "Mulitvariate root finding" << endl;
    wordList multivariateMethods
    (
        rootSolver::dictionaryOneConstructorTablePtr_->toc()
    );
    forAll(multEqns, eqni)
    {
        Info<< "Solving equations: " << nl << multEqns[eqni].name() << endl;
        scalarField x0(2, 0.1);
        const multivariateEquation<scalar>& eqns = multEqns[eqni];
        forAll(multivariateMethods, i)
        {
            dict.set("solver", multivariateMethods[i]);
            autoPtr<rootSolver> solver
            (
                rootSolver::New(eqns, dict)
            );
            Info<< "    roots=" << solver->solve(x0)
                << ", nSteps=" << solver->nSteps()
                << ", error=" << solver->errors() << nl << endl;
        }
    }

    Info<< "done" << endl;

    return 0;
}
