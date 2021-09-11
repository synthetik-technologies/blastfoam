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

    virtual tmp<scalarField> f
    (
        const scalarField& x,
        const label li
    ) const
    {
        tmp<scalarField> fxTmp(new scalarField(x.size()));
        scalarField& fx = fxTmp.ref();
        fx[0] = sqr(x[0]) + sqr(x[1]) - 4.0;
        fx[1] = sqr(x[0]) - x[1] + 1.0;

        return fxTmp;
    }

    virtual void jacobian
    (
        const scalarField& x,
        const label li,
        scalarField& fx,
        scalarSquareMatrix& dfdx
    ) const
    {
        fx = f(x, li);

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

    virtual tmp<scalarField> f
    (
        const scalarField& x,
        const label li
    ) const
    {
        tmp<scalarField> fxTmp(new scalarField(x.size()));
        scalarField& fx = fxTmp.ref();
        fx[0] = x[1] - sqr(x[0]) + x[0];
        fx[1] = sqr(x[0])/16.0 + sqr(x[1]) - 1.0;

        return fxTmp;
    }

    virtual void jacobian
    (
        const scalarField& x,
        const label li,
        scalarField& fx,
        scalarSquareMatrix& dfdx
    ) const
    {
        fx = f(x, li);

        dfdx(0, 0) = -2.0*x[0] + 1.0;
        dfdx(0, 1) = 1.0;
        dfdx(1, 0) = 2.0*x[0]/16.0;
        dfdx(1, 1) = stabilise(2.0*x[1], small);
    }
};


int main(int argc, char *argv[])
{

    PtrList<scalarEquation> uniEqns(2);
    uniEqns.set(0, new testEqn1(0.0, 1.0));
    uniEqns.set(1, new testEqn2(0.0, 1.0));

    PtrList<scalarMultivariateEquation> multEqns(2);
    multEqns.set(0, new multivariateTestEqn1());
    multEqns.set(1, new multivariateTestEqn2());
    dictionary dict;

    Info<< endl;
    Info<< "Univariate root finding" << endl;
    wordList methods(rootSolver::dictionaryConstructorTablePtr_->toc());
    forAll(uniEqns, eqni)
    {
        Info<< "Solving equation " << eqni << endl;

        const scalarEquation& eqn = uniEqns[eqni];
        forAll(methods, i)
        {
            dict.set("solver", methods[i]);
            autoPtr<rootSolver> solver(rootSolver::New(eqn, dict));
            Info<<"roots=" << solver->solve(0.5, 0)
                << ", nSteps=" << solver->nSteps()
                << ", error=" << solver->error() << nl << endl;
        }
    }

    Info<< nl << nl;
    Info<< "Mulitvariate root finding" << endl;
    wordList multivariateMethods
    (
        multivariateRootSolver::dictionaryConstructorTablePtr_->toc()
    );
    forAll(multEqns, eqni)
    {
        Info<< "Solving equation " << eqni << endl;
        scalarField x0(2, 0.0);
        const scalarMultivariateEquation& eqns = multEqns[eqni];
        forAll(multivariateMethods, i)
        {
            dict.set("solver", multivariateMethods[i]);
            autoPtr<multivariateRootSolver> solver
            (
                multivariateRootSolver::New(eqns, dict)
            );
            Info<<"roots=" << solver->solve(x0, 0)
                << ", nSteps=" << solver->nSteps()
                << ", error=" << solver->error() << nl << endl;
        }
    }


    Info<< nl << nl;
    Info<< "Integration of f(x) = cos(x) - x^3 in (0, 0.5)" << endl;
    labelList nSteps({10, 20, 30, 50, 100});
    dict.set("nSteps", 1000);
    const scalar ans = 0.463800538604203;
    wordList intMethods
    (
        scalarIntegrator::dictionaryConstructorTablePtr_->toc()
    );
    forAll(intMethods, i)
    {
        dict.set("integrator", intMethods[i]);
        autoPtr<scalarIntegrator> integrator
        (
            scalarIntegrator::New(uniEqns[0], dict)
        );

        forAll(nSteps, stepi)
        {
            integrator->nSteps() = nSteps[stepi];
            scalar est = integrator->integrate(0, 0.5, 0);
            Info<< "\t" << integrator->nSteps() << " steps: "
                << "estimate=" << est << ", error=" << mag(est - ans) <<endl;
        }
        Info<< endl;
    }

    Info<< nl << "done" << endl;

    return 0;
}
