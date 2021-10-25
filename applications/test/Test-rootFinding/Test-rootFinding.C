#include "dictionary.H"
#include "rootSolver.H"
#include "argList.H"

using namespace Foam;


class testEqn1
:
    public equation
{
public:
    testEqn1(const scalar xMin, const scalar xMax)
    :
        equation(xMin, xMax)
    {}

    virtual ~testEqn1()
    {}

    virtual label nDerivatives() const
    {
        return 2;
    }
    virtual scalar fx(const scalar& x, const label li) const
    {
        return Foam::cos(x) - Foam::pow3(x);
    }
    virtual scalar dfdx(const scalar& x, const label li) const
    {
        return -Foam::sin(x) - 3.0*Foam::sqr(x);
    }
    virtual scalar d2fdx2(const scalar& x, const label li) const
    {
        return -Foam::cos(x) - 6.0*x;
    }
};

class testEqn2
:
    public equation
{
public:
    testEqn2(const scalar xMin, const scalar xMax)
    :
        equation(xMin, xMax)
    {}

    virtual ~testEqn2()
    {}

    virtual label nDerivatives() const
    {
        return 2;
    }
    virtual scalar fx(const scalar& x, const label li) const
    {
        return Foam::exp(x) - 10.0*x;
    }
    virtual scalar dfdx(const scalar& x, const label li) const
    {
        return Foam::exp(x) - 10.0;
    }
    virtual scalar d2fdx2(const scalar& x, const label li) const
    {
        return Foam::exp(x);
    }
};


int main(int argc, char *argv[])
{

    PtrList<equation> uniEqns(2);
    wordList uniEqnStrs(2);
    uniEqns.set(0, new testEqn1(0.0, 1.0));
    uniEqnStrs[0] = "f(x) = cos(x) - x^3";
    uniEqns.set(1, new testEqn2(0.0, 1.0));
    uniEqnStrs[1] = "f(x) = exp(x) - 10*x";

    dictionary dict;

    Info<< endl;
    Info<< "Univariate root finding" << endl;
    wordList methods(rootSolver::dictionaryTwoConstructorTablePtr_->toc());
    forAll(uniEqns, eqni)
    {
        Info<< "Solving " << uniEqnStrs[eqni] << endl;
        const equation& eqn = uniEqns[eqni];
        forAll(methods, i)
        {
            dict.set("solver", methods[i]);
            autoPtr<rootSolver> solver(rootSolver::New(eqn, dict));
            Info<<"    root=" << solver->solve(0.5, 0)
                << ", nSteps=" << solver->nSteps()
                << ", error=" << solver->error() << nl << endl;
        }
    }

    Info<< "done" << endl;

    return 0;
}
