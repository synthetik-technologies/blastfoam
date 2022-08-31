#include "dictionary.H"
#include "linearLeastSquares.H"
#include "nonLinearLeastSquares.H"
#include "CoefficientEquationsFwd.H"
#include "createEquations.H"
#include "UnivariateEquationsFwd.H"
#include "argList.H"
#include "Random.H"

using namespace Foam;

class eqn1
:
    public ScalarCoefficientEquation
{
public:

    eqn1()
    :
        ScalarCoefficientEquation(3)
    {}

    label nDerivatives() const
    {
        return 0;
    }

    scalar fx(const scalar x, const label li) const
    {
        return
            coeffs_[0]
          + coeffs_[1]*x
          + coeffs_[2]*sqr(x);
    }

    virtual void coeffJ
    (
        const XType& x,
        const label i,
        RectangularMatrix<scalar>& J,
        const label li
    ) const
    {
        J(i, 0) = -1.0;
        J(i, 1) = -x[0];
        J(i, 2) = -sqr(x[0]);
    }
};

class eqn2
:
    public ScalarUnivariateCoefficientEquation
{
public:

    eqn2()
    :
        ScalarUnivariateCoefficientEquation(2, 6)
    {}

    label nDerivatives() const 
    {
        return 0;
    }

    scalar fX(const VarType& x, const label li) const
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
        const XType& x,
        const label i,
        RectangularMatrix<scalar>& J,
        const label li
    ) const
    {
        J(i, 0) = -1.0;
        J(i, 1) = -x[0];
        J(i, 2) = -x[1];
        J(i, 3) = -sqr(x[0]);
        J(i, 4) = -x[0]*x[1];
        J(i, 5) = -sqr(x[1]);
    }
};

int main(int argc, char *argv[])
{
    Random rand(0);
    const label n = 100;
    const scalar xMax = 20.0;
    List<scalar> x1(n);
    List<vector2D> x2(n);
    scalarField y1(n);
    scalarField y2(n);
    scalarField y3(n);
    scalarField y4(n);
    forAll(x1, i)
    {
        x1[i] = rand.scalarAB(0, xMax);
        y1[i] = 0.5 + x1[i]*1.7;

        x2[i][0] = x1[i];
        x2[i][1] = rand.scalarAB(0, xMax);
        y2[i] =
            y1[i]
          - x2[i][1]*3.5;

        y3[i] = y1[i] - 0.5*sqr(x1[i]);

        y4[i] =
            y2[i]
          + sqr(x2[i][0])*1.5
          + x2[i][0]*x2[i][1]*4.2
          - sqr(x2[i][1])*3.5;
    }

    dictionary dict;
    Info<< nl << "leastSquares:" << endl;

    {
        ScalarLinearEquation eqn;
        linearLeastSquares solver;
        solver.findCoeffs(eqn, x1, y1);
        Info<<eqn.coeffs()<<endl;
    }
    {
        linearLeastSquares solver;
        autoPtr<scalarUnivariateEquation> eqn(solver.createEquation(x2, y2));
        Info<<dynamic_cast<const scalarCoefficients&>(eqn()).coeffs()<<endl;
    }
    {
        eqn1 eqn;
        nonLinearLeastSquares solver;
        solver.findCoeffs(eqn, x2, y3);
        Info<<eqn.coeffs()<<endl;
    }
    {
        eqn2 eqn;
        nonLinearLeastSquares solver;
        solver.findCoeffs(eqn, x2, y4);
        Info<<eqn.coeffs()<<endl;
    }
    Info<<nl;

    Info<< nl << "done" << endl;

    return 0;
}
