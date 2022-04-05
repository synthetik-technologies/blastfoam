#include "dictionary.H"
#include "minimizationScheme.H"
#include "EquationsFwd.H"
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
            "f(x, y, z) = ((x + 2)*sin(x))^2 + (y + 2)^2 + (z + 3)^2",
            scalarField(scalarList{2, -4, -4}),
            scalarField(scalarList{4, 4, 4})
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
        return
            sqr((x[0] + 2.0)*Foam::sin(x[0])) + sqr(x[1] + 2.0) + sqr(x[2] + 3.0);
    }
};

class testEqn2
:
    public ScalarEquation
{
public:
    testEqn2()
    :
        ScalarEquation("f(x) = (x - 2)^2", -3, 3)
    {}

    virtual ~testEqn2()
    {}

    virtual label nDerivatives() const
    {
        return 0;
    }
    virtual scalar fx(const scalar x, const label li) const
    {
        return sqr(x - 2.0);
    }
    virtual scalar dfdx(const scalar x, const label li) const
    {
        return 2.0*(x - 2.0);
    }
};


int main(int argc, char *argv[])
{
    PtrList<scalarUnivariateEquation> multEqns(2);
    multEqns.set(0, new testEqn1());
    multEqns.set(1, new testEqn2());

    dictionary dict;
//     dict.add("maxSteps", 10);
    dict.add("nParticles", 1000);

    Info<< "Minimization" << endl;
    minimizationScheme::debug = 0;
    forAll(multEqns, eqni)
    {
        Info<< nl << "Solving equations: " << nl << multEqns[eqni].name() << endl;
        scalarField x0Orig(scalarList({0, 0, 0}));
        const scalarUnivariateEquation& eqns = multEqns[eqni];

        wordList multivariateMethods
        (
            eqns.nVar() == 1
          ? minimizationScheme::dictionaryUnivariateConstructorTablePtr_->toc()
          : minimizationScheme::dictionaryMultivariateConstructorTablePtr_->toc()
        );

        forAll(multivariateMethods, i)
        {
            scalarField x0(x0Orig);
            x0.resize(eqns.nVar());
            dict.set("solver", multivariateMethods[i]);
            autoPtr<minimizationScheme> solver
            (
                minimizationScheme::New(eqns, dict)
            );
            scalarField localMin(solver->solve());
            Info<< "    " << multivariateMethods[i] << ": "
                << "Local minimums = " << localMin << endl;
        }
    }

    Info<< nl << "done" << endl;

    return 0;
}
