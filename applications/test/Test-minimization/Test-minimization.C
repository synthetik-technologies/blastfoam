#include "dictionary.H"
#include "minimizationScheme.H"
#include "equation.H"
#include "argList.H"

using namespace Foam;


class testEqn1
:
    public ScalarEquation<scalarField>
{
public:
    // Constructors
    testEqn1()
    :
        ScalarEquation<scalarField>
        (
            3,
            scalarField(3, -30.0),
            scalarField(3, 30.0)
        )
    {}

    //- Destructor
    virtual ~testEqn1()
    {}

    virtual label nDerivatives() const
    {
        return 0;
    }
    virtual void f
    (
        const scalarField& x,
        const label li,
        scalar& fx
    ) const
    {
        fx = sqr((x[0] + 2.0)*Foam::sin(x[0])) + sqr(x[1] + 2.0) + sqr(x[2] + 3.0);
        // fx[0] = x[0]*sqr(x[1]);
    }
};

class testEqn2
:
    public equation
{
public:
    testEqn2()
    :
        equation(-10, 10)
    {}

    virtual ~testEqn2()
    {}

    virtual label nDerivatives() const
    {
        return 0;
    }
    virtual scalar fx(const scalar& x, const label li) const
    {
        return sqr(x - 2.0);
    }
};


int main(int argc, char *argv[])
{
    PtrList<scalarEquation> multEqns(2);
    wordList multiEqnStrs(2);
    multEqns.set(0, new testEqn1());
    multiEqnStrs[0] =
        word("    f(x1, x2) = ((x1 + 1)*sin(x))^2 + (x2 + 2)^2 + (x3 + 3)^3\n");
    multEqns.set(1, new testEqn2());
    multiEqnStrs[1] =
        word("    f(x1) = (x1 - 2.0^2)\n");

    dictionary dict;
    // dict.add("maxSteps", 10);

    Info<< "Minimization" << endl;
    wordList multivariateMethods
    (
        minimizationScheme::dictionaryOneConstructorTablePtr_->toc()
    );
    minimizationScheme::debug = 1;
    forAll(multEqns, eqni)
    {
        Info<< nl << "Solving equations: " << nl << multiEqnStrs[eqni] << endl;
        scalarField x0Orig(scalarList({2, 2, 2}));
        const scalarEquation& eqns = multEqns[eqni];
        forAll(multivariateMethods, i)
        {
            scalarField x0(x0Orig);
            x0.resize(eqns.nVar());
            dict.set("solver", multivariateMethods[i]);
            autoPtr<minimizationScheme> solver
            (
                minimizationScheme::New(eqns, dict)
            );
            scalarField localMin(solver->solve(x0));
            Info<< "    Local minimums=" << localMin << endl;
        }
    }

    Info<< nl << "done" << endl;

    return 0;
}
