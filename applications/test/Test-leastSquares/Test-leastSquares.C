#include "dictionary.H"
#include "linearLeastSquares.H"
#include "nonLinearLeastSquares.H"
#include "createEquations.H"
#include "UnivariateEquationsFwd.H"
#include "argList.H"
#include "Random.H"

using namespace Foam;

class eqn1
:
    public ScalarUnivariateEquation,
    public coefficients<scalar, scalar>
{
public:

    eqn1()
    :
        ScalarUnivariateEquation(2, {-great, -great}, {great, great}),
        coefficients<scalar, scalar>(3)
    {}

    label nDerivatives() const 
    {
        return 0;
    }

    scalar fX(const scalarList& x, const label li) const
    {
        scalar res = coeffs_[0];
        forAll(x, i)
        {
            res += coeffs_[i+1]*x[i];
        }
        return res;
    }
};
class eqn2
:
    public ScalarUnivariateEquation,
    public coefficients<scalar, scalar>
{
public:

    eqn2()
    :
        ScalarUnivariateEquation(2, {-great, -great}, {great, great}),
        coefficients<scalar, scalar>(6)
    {}

    label nDerivatives() const 
    {
        return 0;
    }

    scalar fX(const scalarList& x, const label li) const
    {
        return 
            coeffs_[0]
          + coeffs_[1]*x[0] 
          + coeffs_[2]*x[1]
          + coeffs_[3]*sqr(x[0]) 
          + coeffs_[4]*x[0]*x[1]
          + coeffs_[5]*sqr(x[1]);
    }

    virtual void coeffJ
    (
        const List<scalarList>& x,
        const label li,
        RectangularMatrix<scalar>& J
    ) const
    {
        J.setSize(x.size(), 6);
        forAll(x, i)
        {
            J(i, 0) = -1.0;
            J(i, 1) = -x[i][0];
            J(i, 2) = -x[i][1];
            J(i, 3) = -sqr(x[i][0]);
            J(i, 4) = -x[i][0]*x[i][1];
            J(i, 5) = -sqr(x[i][1]);
        }
    }
};

int main(int argc, char *argv[])
{
    Random rand(0);
    const label n = 100;
    const scalar xMax = 20.0;
    List<scalarList> x(n, scalarList(2));
    scalarList y(n);
    forAll(x, i)
    {
        x[i][0] = rand.scalarAB(0, xMax);
        x[i][1] = rand.scalarAB(0, xMax);
        y[i] = 
            x[i][1]*(x[i][0])*1.7
          - sqr(x[i][1] - 5.0)
          + rand.scalarAB(-xMax, xMax);
    }

    dictionary dict;
    Info<< nl << "leastSquares:" << endl;

    {
        eqn1 eqn;
        linearLeastSquares solver;
        solver.findCoeffs(eqn, x, y);
        Info<<eqn.coeffs()<<endl;
    }
    Info<<nl;

    {
        eqn2 eqn;
        nonLinearLeastSquares solver;
        solver.findCoeffs(eqn, x, y);
        Info<<eqn.coeffs()<<endl;
    }
    Info<<nl;

    Info<< nl << "done" << endl;

    return 0;
}
