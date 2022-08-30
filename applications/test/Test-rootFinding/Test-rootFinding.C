#include "dictionary.H"
#include "univariateRootSolver.H"
#include "argList.H"
#include "createEquations.H"

using namespace Foam;


createNamedEquation2
(
    testEqn1,
    "f(x) = cos(x) - x^3",
    scalar,
    0.01, 1.0,
    Foam::cos(x) - Foam::pow3(x),
    -Foam::sin(x) - 3.0*Foam::sqr(x),
    -Foam::cos(x) - 6.0*x
);

createNamedEquation2
(
    testEqn2,
    "f(x) = exp(x) - 10*x",
    scalar,
    0.01, 1.0,
    Foam::exp(x) - 10.0*x,
    Foam::exp(x) - 10.0,
    Foam::exp(x)
);

class testEqn3
:
    public ScalarMultivariateEquation
{
public:
    // Constructors
    testEqn3()
    :
        ScalarMultivariateEquation
        (
            2,
            {0.0, 0.0},
            {2.0, 2.0},
            {
                "f1(x1, x2) = x1^2 + x2^2 - 4.0",
                "f2(x1, x2) = x1^2 - x2 + 1.0"
            }
        )
    {}

    //- Destructor
    virtual ~testEqn3()
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
        const UList<scalar>& x,
        const label li,
        scalarList& fx
    ) const
    {
        fx[0] = sqr(x[0]) + sqr(x[1]) - 4.0;
        fx[1] = sqr(x[0]) - x[1] + 1.0;
    }
    virtual void jacobian
    (
        const UList<scalar>& x,
        const label li,
        scalarList& fx,
        RectangularMatrix<scalar>& dfdx
    ) const
    {
        FX(x, li, fx);
        dfdx(0, 0) = stabilise(2.0*x[0], small);
        dfdx(0, 1) = 2.0*x[1];
        dfdx(1, 0) = 2.0*x[0];
        dfdx(1, 1) = -1.0;
    }
};


class testEqn4
:
    public ScalarMultivariateEquation
{
public:
    // Constructors
    testEqn4()
    :
        ScalarMultivariateEquation
        (
            2,
            {-10.0, -10.0},
            {10.0, 10.0},
            {
                "f1(x1, x2) = x2^2 + x1^2 + x1",
                "f2(x1, x2) = (x1^2)/16 - x2^2 - 1.0"
            }
        )
    {}

    //- Destructor
    virtual ~testEqn4()
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
        const UList<scalar>& x,
        const label li,
        scalarList& fx
    ) const
    {
        fx[0] = x[1] - sqr(x[0]) + x[0];
        fx[1] = sqr(x[0])/16.0 + sqr(x[1]) - 1.0;
    }
    virtual void jacobian
    (
        const UList<scalar>& x,
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

    PtrList<scalarMultivariateEquation> eqns(4);
    eqns.set(0, new testEqn1());
    eqns.set(1, new testEqn2());
    eqns.set(2, new testEqn3());
    eqns.set(3, new testEqn4());

    dictionary dict;

    Info<< "*****************************************" << nl
        << "Testing root finders" << nl
        << "*****************************************" << nl << endl;
    List<wordList> methods =
    {
        univariateRootSolver::dictionaryTwoConstructorTablePtr_->toc(),
        univariateRootSolver::dictionaryTwoConstructorTablePtr_->toc(),
        rootSolver::dictionaryOneConstructorTablePtr_->toc(),
        rootSolver::dictionaryOneConstructorTablePtr_->toc()
    };
    List<scalarField> x0 =
    {
        scalarField(1, 0.5),
        scalarField(1, 0.5),
        scalarField(2, 0.1),
        scalarField(scalarList({-0.1, 0.1}))
    };

    List<scalarList> roots =
    {
        {0.865474033101615},
        {0.111832559158963},
        {0.889543928069176, 1.791287856740247},
        {-0.612743189209448, 0.988197405130338}
    };

    rootSolver::debug = 2;
    forAll(eqns, eqni)
    {
        Info<< "*****************************************" << nl
            << "Solving:" << nl << eqns[eqni].eqnString() << nl;
        if (roots[eqni].size() > 1)
        {
            Info<< "roots=" << roots[eqni] <<nl;
        }
        else
        {
            Info<< "root=" << roots[eqni][0] <<nl;
        }
        Info<< "*****************************************" << endl;
        const scalarMultivariateEquation& eqn = eqns[eqni];
        forAll(methods[eqni], i)
        {
            autoPtr<rootSolver> solver
            (
                rootSolver::New(methods[eqni][i], eqn, dict)
            );
            incrIndent(Info);
            scalarField r(solver->solve(x0[eqni], 0));
            if (roots[eqni].size() > 1)
            {
                Info<< indent << "error=" << mag(r-roots[eqni]) <<nl;
            }
            else
            {
                Info<< indent<< "error=" << mag(r[0]-roots[eqni][0]) <<nl;
            }
            Info << decrIndent << endl;
        }
    }

    Info<< "done" << endl;

    return 0;
}
