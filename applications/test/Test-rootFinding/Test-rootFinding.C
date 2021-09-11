#include "dictionary.H"
#include "rootSolver.H"
#include "multivariateRootSolver.H"
#include "argList.H"
#include "IntegratorsFwd.H"

using namespace Foam;


class testEqn
:
    public scalarEquation
{
public:
    // Constructor
    testEqn()
    {}

    testEqn(const scalar xMin, const scalar xMax)
    :
        scalarEquation(xMin, xMax)
    {}

    //- Destructor
    virtual ~testEqn()
    {}


    // Member Functions

        //- Return the function value
        virtual scalar f(const scalar x, const label li) const
        {
            return Foam::cos(x) - Foam::pow3(x);
        }

        //- Calculate the first derivative of the equation
        virtual scalar dfdx(const scalar x, const label li) const
        {
            return -Foam::sin(x) - 3.0*Foam::sqr(x);
        }

        //- Calculate the second derivative of the equation
        virtual scalar d2fdx2(const scalar x, const label li) const
        {
            return -Foam::cos(x) - 6.0*x;
        }
};


class multivariateTestEqn
:
    public scalarMultivariateEquation
{
public:
    // Constructors
    multivariateTestEqn()
    {}

    multivariateTestEqn(const scalarField& xMins, const scalarField& xMaxs)
    :
        scalarMultivariateEquation(xMins, xMaxs)
    {}

    //- Destructor
    virtual ~multivariateTestEqn()
    {}


    // Member Functions

        //- Return the number of equations in the system
        virtual label nEqns() const
        {
            return 2;
        }

        //- Return the function value
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

        //- Calculate the first derivative of the equation
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

int main(int argc, char *argv[])
{

    testEqn eqn(0, 2.0);
    multivariateTestEqn eqns;
    dictionary dict;

    Info<< endl;
    Info<< "Univariate root finding" << endl;
    wordList methods(rootSolver::dictionaryConstructorTablePtr_->toc());
    forAll(methods, i)
    {
        dict.set("solver", methods[i]);
        autoPtr<rootSolver> solver(rootSolver::New(eqn, dict));
        Info<<"roots=" << solver->solve(0.5, 0)
            << ", nSteps=" << solver->nSteps()
            << ", error=" << solver->error() << nl << endl;
    }

    Info<< nl << nl;
    Info<< "Mulitvariate root finding" << endl;
    wordList multivariateMethods
    (
        multivariateRootSolver::dictionaryConstructorTablePtr_->toc()
    );
    forAll(multivariateMethods, i)
    {
        dict.set("solver", multivariateMethods[i]);
        autoPtr<multivariateRootSolver> solver
        (
            multivariateRootSolver::New(eqns, dict)
        );
        scalarField x0(2, 0.0);
        Info<<"root=" << solver->solve(x0, 0)
            << ", nSteps=" << solver->nSteps()
            << ", errors=" << solver->errors() << nl << endl;
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
            scalarIntegrator::New(eqn, dict)
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
