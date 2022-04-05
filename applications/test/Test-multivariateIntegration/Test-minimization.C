#include "dictionary.H"
#include "MultivariateIntegratorsFwd.H"
#include "UnivariateEquationsFwd.H"
#include "argList.H"

using namespace Foam;


class testEqn1
:
    public ScalarUnivariateEquation
{
public:
    // Constructors
    testEqn1()
    :
        ScalarUnivariateEquation
        (
            "f(x, y, z) = x*y*z",
            {-great, -great, -great},
            {great, great, great}
        )
    {}

    //- Destructor
    virtual ~testEqn1()
    {}

    virtual label nDerivatives() const
    {
        return 0;
    }
    virtual scalar fX
    (
        const scalarList& x,
        const label li
    ) const
    {
        return sqr(x[0] + Foam::exp(x[1]))*x[1]*x[2] - Foam::cos(x[2]);
    }
};

class testEqn2
:
    public ScalarUnivariateEquation
{
public:
    // Constructors
    testEqn2()
    :
        ScalarUnivariateEquation
        (
            "f(x) = (x - 2.0^2)",
            {-great},
            {great}
        )
    {}

    //- Destructor
    virtual ~testEqn2()
    {}

    virtual label nDerivatives() const
    {
        return 0;
    }
    virtual scalar fX
    (
        const scalarList& x,
        const label li
    ) const
    {
        return x[0] - sqr(2.0);
    }
};

int main(int argc, char *argv[])
{
    PtrList<scalarUnivariateEquation> multEqns(1);
    multEqns.set(0, new testEqn1());
    multEqns.set(1, new testEqn2());

    dictionary dict;

    Info<< "Multivariate integration" << endl;
    wordList intMethods1
    (
        scalarMultivariateIntegrator::dictionaryConstructorTablePtr_->toc()
    );

    Info<< "Integration of " << multEqns[0].name() << " from (0, 0, 0) : (1, 1, 1)"
        << endl;

    scalarField x0(3, 0.0);
    scalarField x1(3, 1.0);
//     const scalar ans1 = 0.125;
//     const scalar ans1 = 8.333333333333334e-02;
    const scalar ans1 = 0.790494360891768;
    forAll(intMethods1, i)
    {
        dict.set("integrator", intMethods1[i]);
        dict.set("tolerance", scalarField(3, 1e-6));
        dict.set("nNodes", labelList(3, 10));
//         dict.set("adaptive", "false");
        dict.set("maxSplits", labelList(3, 1));
        autoPtr<scalarMultivariateIntegrator> integrator
        (
            scalarMultivariateIntegrator::New(multEqns[0], dict)
        );
        scalar est = integrator->integrate(x0, x1, 0);
        Info<< "\t" << integrator->nEvals() << " adaptive intervals: "
            << "estimate=" << est << ", error=" << mag(est - ans1) <<endl;
//         forAll(nIntervals, inti)
//         {
//             integrator->setNIntervals(nIntervals[inti]);
//             est = integrator->integrate(0.0, 2.5, 0);
//             Info<< "\t" << integrator->nIntervals() << " steps: "
//                 << "estimate=" << est << ", error=" << mag(est - ans1) <<endl;
//         }
        Info<< endl;
    }

    Info<< nl << "done" << endl;

    return 0;
}
