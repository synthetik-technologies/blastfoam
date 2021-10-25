#include "dictionary.H"
#include "argList.H"
#include "UnivariateEquation.H"
#include "IntegratorsFwd.H"

using namespace Foam;


class testEqn1
:
    public UnivariateEquation<scalar>
{
public:
    testEqn1(const scalar xMin, const scalar xMax)
    :
        UnivariateEquation<scalar>(1, xMin, xMax)
    {}

    virtual ~testEqn1()
    {}

    virtual label nDerivatives() const
    {
        return 2;
    }

    virtual scalar fx(const scalar x, const label li) const
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
    public UnivariateEquation<scalar>
{
public:
    testEqn2(const scalar xMin, const scalar xMax)
    :
        UnivariateEquation<scalar>(1, xMin, xMax)
    {}

    virtual ~testEqn2()
    {}

    virtual label nDerivatives() const
    {
        return 2;
    }

    virtual scalar fx(const scalar x, const label li) const
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

int main(int argc, char *argv[])
{

    testEqn1 eqn(0.0, 1.0);
    dictionary dict;

    Info<< "Integration of f(x) = cos(x) - x^3 in (0, 2.5)" << endl;
    labelList nIntervals({1, 5, 10, 25, 50});
    const scalar ans = -9.167152855896045;
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

        forAll(nIntervals, inti)
        {
            integrator->setNIntervals(nIntervals[inti]);
            scalar est = integrator->integrate(0, 2.5, 0);
            Info<< "\t" << integrator->nIntervals() << " steps: "
                << "estimate=" << est << ", error=" << mag(est - ans) <<endl;
        }
        Info<< endl;
    }

    Info<< nl << "done" << endl;

    return 0;
}
