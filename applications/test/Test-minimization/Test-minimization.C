#include "dictionary.H"
#include "minimizationScheme.H"
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

    virtual label nDerivatives() const
    {
        return 2;
    }

    virtual scalar f(const scalar x, const label li) const
    {
        return mag(x - 2.0) + sqr(x - 1.0);
    }
    virtual scalar dfdx(const scalar x, const label li) const
    {
        return (x - 2.0)/max(mag(x - 2.0), small) + 2.0*(x - 1.0);
    }
    virtual scalar d2fdx2(const scalar x, const label li) const
    {
        return 2.0;
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

    virtual label nDerivatives() const
    {
        return 2;
    }

    virtual scalar f(const scalar x, const label li) const
    {
        return sqr(x - 2.0);
    }
    virtual scalar dfdx(const scalar x, const label li) const
    {
        return 2.0*(x - 2.0);
    }
    virtual scalar d2fdx2(const scalar x, const label li) const
    {
        return 2.0;
    }
};


int main(int argc, char *argv[])
{

    PtrList<scalarEquation> uniEqns(2);
    wordList uniEqnStrs(2);
    uniEqns.set(0, new testEqn1(1.0, 3.0));
    uniEqnStrs[0] = "f(x) = cos(x) - x^3";
    uniEqns.set(1, new testEqn2(1.0, 5.0));
    uniEqnStrs[1] = "f(x) = exp(x) - 10*x";

    dictionary dict;

    Info<< endl;
    Info<< "Minimization" << endl;
    wordList methods
    (
        minimizationScheme::dictionaryTwoConstructorTablePtr_->toc()
    );
    dict.set("tolerance", 1e-5);
    forAll(uniEqns, eqni)
    {
        Info<< "Solving " << uniEqnStrs[eqni] << endl;
        incrIndent(Info);

        const scalarEquation& eqn = uniEqns[eqni];
        forAll(methods, i)
        {
            dict.set("solver", methods[i]);
            autoPtr<minimizationScheme> solver
            (
                minimizationScheme::New(eqn, dict)
            );
            Info<<"\tmin=" << solver->solve()
                << ", nSteps=" << solver->nSteps()
                << ", error=" << solver->error() << nl << endl;
        }
        Info<<nl;
        decrIndent(Info);
    }

    Info<< nl << "done" << endl;

    return 0;
}
