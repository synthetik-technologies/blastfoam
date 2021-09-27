#include "dictionary.H"
#include "rootSolver.H"
#include "multivariateRootSolver.H"
#include "argList.H"
#include "IntegratorsFwd.H"

using namespace Foam;


class testEqn1
:
    public scalarEquation
{
public:
    testEqn1(const scalar xMin, const scalar xMax)
    :
        scalarEquation(xMin, xMax)
    {}

    virtual ~testEqn1()
    {}

    virtual scalar f(const scalar x, const label li) const
    {
        return Foam::cos(x) - Foam::pow3(x);
    }
    virtual scalar dfdx(const scalar x, const label li) const
    {
        return -Foam::sin(x) - 3.0*Foam::sqr(x);
    }
    virtual scalar d2fdx2(const scalar x, const label li) const
    {
        return -Foam::cos(x) - 6.0*x;
    }
};

class testEqn2
:
    public scalarEquation
{
public:
    testEqn2(const scalar xMin, const scalar xMax)
    :
        scalarEquation(xMin, xMax)
    {}

    virtual ~testEqn2()
    {}

    virtual scalar f(const scalar x, const label li) const
    {
        return Foam::exp(x) - 10.0*x;
    }
    virtual scalar dfdx(const scalar x, const label li) const
    {
        return Foam::exp(x) - 10.0;
    }
    virtual scalar d2fdx2(const scalar x, const label li) const
    {
        return Foam::exp(x);
    }
};


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

    PtrList<scalarEquation> uniEqns(2);
    wordList uniEqnStrs(2);
    uniEqns.set(0, new testEqn1(0.0, 1.0));
    uniEqnStrs[0] = "f(x) = cos(x) - x^3";
    uniEqns.set(1, new testEqn2(0.0, 1.0));
    uniEqnStrs[1] = "f(x) = exp(x) - 10*x";

    PtrList<scalarMultivariateEquation> multEqns(2);
    wordList multiEqnStrs(2);
    multEqns.set(0, new multivariateTestEqn1());
    multiEqnStrs[0] =
        word("\tf1(x1, x2) = x1^2 + x2^2 - 4.0\n")
      + word("\tf2(x1, x2) = x1^2 - x2 + 1.0\n");
    multEqns.set(1, new multivariateTestEqn2());
    multiEqnStrs[1] =
        word("\tf1(x1, x2) = x2^2 + x1^2 + x1\n")
      + word("\tf2(x1, x2) = (x1^2)/16 - x2^2 - 1.0\n");

    dictionary dict;

    Info<< endl;
    Info<< "Univariate root finding" << endl;
    wordList methods(rootSolver::dictionaryConstructorTablePtr_->toc());
    forAll(uniEqns, eqni)
    {
        Info<< "Solving " << uniEqnStrs[eqni] << endl;
        incrIndent(Info);

        const scalarEquation& eqn = uniEqns[eqni];
        forAll(methods, i)
        {
            dict.set("solver", methods[i]);
            autoPtr<rootSolver> solver(rootSolver::New(eqn, dict));
            Info<<"\troot=" << solver->solve(0.5, 0)
                << ", nSteps=" << solver->nSteps()
                << ", error=" << solver->error() << nl << endl;
        }
        Info<<nl;
        decrIndent(Info);
    }

    Info<< nl << nl;
    Info<< "Mulitvariate root finding" << endl;
    wordList multivariateMethods
    (
        multivariateRootSolver::dictionaryConstructorTablePtr_->toc()
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
            Info<<"\troots=" << solver->solve(x0, 0)
                << ", nSteps=" << solver->nSteps()
                << ", error=" << solver->error() << nl << endl;
        }
        Info<<nl;
    }


//     Info<< nl << nl;
//     Info<< "Integration of f(x) = cos(x) - x^3 in (0, 2.5)" << endl;
//     labelList nIntervals({1, 5, 10, 25, 50});
//     const scalar ans = -9.167152855896045;
//     wordList intMethods
//     (
//         scalarIntegrator::dictionaryConstructorTablePtr_->toc()
//     );
//     forAll(intMethods, i)
//     {
//         dict.set("integrator", intMethods[i]);
//         autoPtr<scalarIntegrator> integrator
//         (
//             scalarIntegrator::New(uniEqns[0], dict)
//         );
//
//         forAll(nIntervals, inti)
//         {
//             integrator->setNIntervals(nIntervals[inti]);
//             scalar est = integrator->integrate(0, 2.5, 0);
//             Info<< "\t" << integrator->nIntervals() << " steps: "
//                 << "estimate=" << est << ", error=" << mag(est - ans) <<endl;
//         }
//         Info<< endl;
//     }

    Info<< nl << "done" << endl;

    return 0;
}
