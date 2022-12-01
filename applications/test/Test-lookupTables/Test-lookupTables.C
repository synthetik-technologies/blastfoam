#include "dictionary.H"
#include "lookupTables1D.H"
#include "lookupTables2D.H"
#include "lookupTables3D.H"
#include "OFstream.H"
#include "IFstream.H"
#include "Random.H"

#include "univariateRootSolver.H"
#include "EquationsFwd.H"

#include "List2D.H"
#include "List3D.H"

#include "argList.H"
#include "IOmanip.H"
#include "PtrList2D.H"

using namespace Foam;

scalar func1(const scalar x)
{
    return 3.0*x*x + 1.0;
}
scalar dfunc1dx(const scalar x)
{
    return 6.0*x;
}
scalar d2func1dx2(const scalar x)
{
    return 6.0;
}

scalar func2(const scalar x, const scalar y)
{
    return 3.0*x*x + sqr(y)*x;
}
scalar dfunc2dx(const scalar x, const scalar y)
{
    return 6.0*x + sqr(y);
}
scalar dfunc2dy(const scalar x, const scalar y)
{
    return 2.0*y*x;
}
scalar d2func2dx2(const scalar x, const scalar y)
{
    return 6.0;
}
scalar d2func2dy2(const scalar x, const scalar y)
{
    return 2.0*x;
}
scalar d2func2dxdy(const scalar x, const scalar y)
{
    return 2.0*y;
}


scalar func3(const scalar x, const scalar y, const scalar z)
{
    return 3.0*x*x*y*z + y*y*Foam::sqr(z);
}
scalar dfunc3dx(const scalar x, const scalar y, const scalar z)
{
    return 6.0*x*y*z;
}
scalar dfunc3dy(const scalar x, const scalar y, const scalar z)
{
    return 3.0*x*x*z + 2.0*y*Foam::sqr(z);
}
scalar dfunc3dz(const scalar x, const scalar y, const scalar z)
{
    return 3.0*x*x*y + 2.0*y*y*z;
}
scalar d2func3dx2(const scalar x, const scalar y, const scalar z)
{
    return 6*y*z;
}
scalar d2func3dy2(const scalar x, const scalar y, const scalar z)
{
    return 2.0*sqr(z);
}
scalar d2func3dz2(const scalar x, const scalar y, const scalar z)
{
    return 2.0*y*y;
}
scalar d2func3dxdy(const scalar x, const scalar y, const scalar z)
{
    return 6.0*x*z;
}
scalar d2func3dxdz(const scalar x, const scalar y, const scalar z)
{
    return 6.0*x*y;
}
scalar d2func3dydz(const scalar x, const scalar y, const scalar z)
{
    return 3.0*x*x + 4.0*y*z;
}

template<class Type>
void print(const word& name, const Type& x, const Type& ans)
{
    Info<< name << " (calc/true): " << x <<  "/" << Foam::name(ans)
        << ", error (abs, rel): " << mag(x - ans) << '/' << mag(x - ans)/max(mag(ans), small)
        << endl;
}
int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    IFstream is("tableDict");
    dictionary dict(is);

    vectorLookupTable3D t3;
    vectorLookupTable2D t2;
    // Create some tables
    label nx = dict.subDict("table1D").lookup<label>("nX");
    label ny = 35;
    label nz = 25;
    scalar xMin = dict.subDict("table1D").lookup<scalar>("minX");
    scalar yMin = 0.1;
    scalar zMin = 0.001;
    scalar xMax = dict.subDict("table1D").lookup<scalar>("maxX");
    scalar yMax = 3.0;
    scalar zMax = 4.0;
    scalar dx = (xMax - xMin)/scalar(nx);
    scalar dy = (yMax - yMin)/scalar(ny);
    scalar dz = (zMax - zMin)/scalar(nz);

    scalarField x(nx+1);
    scalarField y(ny+1);
    scalarField z(nz+1);

    {
        OFstream outX("x.csv");
        for (label i = 0; i <= nx; i++)
        {
            x[i] =  xMin + dx*scalar(i);
            outX << x[i] << ";";
        }
    }
    {
        OFstream outY("y.csv");
        for (label j = 0; j <= ny; j++)
        {
            y[j] =  yMin + dy*scalar(j);
            outY << y[j] << ";\n";
        }
    }
    {
        OFstream outZ("z.csv");
        outZ << "# abc"<<endl;
        for (label k = 0; k <= nz; k++)
        {
            z[k] =  zMin + dz*scalar(k);
            outZ << z[k] << "\n";
        }
    }

    {
        OFstream out1D("table1D.csv");
        for (label i = 0; i <= nx; i++)
        {
            out1D << func1(x[i]);
            if ( i != nx)
            {
                out1D << ",";
            }
        }
    }
    {
        OFstream out2D("table2D.csv");
        for (label i = 0; i <= nx; i++)
        {
            for (label j = 0; j <= ny; j++)
            {
                out2D << func2(x[i], y[j]);
                if ( j != ny)
                {
                    out2D << ",";
                }
            }
            out2D <<endl;
        }
    }
    {
        OFstream out3D("table3D.csv");
        for (label k = 0; k <= nz; k++)
        {
            for (label j = 0; j <= ny; j++)
            {
                for (label i = 0; i <= nx; i++)
                {
                    out3D << func3(x[i], y[j], z[k]);
                    if ( i != nx)
                    {
                        out3D << ",";
                    }
                }
                out3D <<";";
            }
            out3D<<endl;
        }
    }

    {
        Random rand(0);
        label nSamples = 10000;

        OFstream outX("sparseX.csv");
        OFstream outY("sparseY.csv");
        OFstream outZ("sparseZ.csv");
        OFstream outF2("sparseF2.csv");
        OFstream outF3("sparseF3.csv");
        for (label i = 0; i < nSamples; i++)
        {
            scalar x = rand.scalarAB(xMin, xMax);
            scalar y = rand.scalarAB(yMin, yMax);
            scalar z =  rand.scalarAB(zMin, zMax);
            outX << x << endl;
            outY << y << endl;
            outZ << z << endl;
            outF2 << func2(x, y) << endl;
            outF3 << func3(x, y, z) << endl;
        }
    }



    scalar xTest = 1.435;
    scalar yTest = 1.3346;
    scalar zTest = 2.5676;

    Info<<nl<<"1D table:" << endl;
    scalarLookupTable1D table1(dict.subDict("table1D"), "x", "f");
    scalar xFound = table1.reverseLookup(table1.lookup(xTest));
    print("f", table1.lookup(xTest), func1(xTest));
    print("dfdx", table1.dFdX(xTest), dfunc1dx(xTest));
    print("d2fdx2", table1.d2FdX2(xTest), d2func1dx2(xTest));
    print("reverse", xFound, xTest);

    Info<<nl<<"2D table:" << endl;
    lookupTable2D<scalar> table2(dict.subDict("table2D"), "x", "y", "f");
    xFound = table2.reverseLookupX(table2.lookup(xTest, yTest), yTest);
    scalar yFound = table2.reverseLookupY(table2.lookup(xTest, yTest), xTest);
    print("f", table2.lookup(xTest, yTest), func2(xTest, yTest));
    print("dfdx", table2.dFdX(xTest, yTest), dfunc2dx(xTest, yTest));
    print("dfdy", table2.dFdY(xTest, yTest), dfunc2dy(xTest, yTest));
    print("d2fdx2", table2.d2FdX2(xTest, yTest), d2func2dx2(xTest, yTest));
    print("d2fdy2", table2.d2FdY2(xTest, yTest), d2func2dy2(xTest, yTest));
    print("d2fdxdy", table2.d2FdXdY(xTest, yTest), d2func2dxdy(xTest, yTest));
    print("reverseX", xFound, xTest);
    print("reverseY", yFound, yTest);

    Info<<nl<<"2D table from least squares:" << endl;
    lookupTable2D<scalar> table2_ls(dict.subDict("table2D_ls"), "x", "y", "f");
    xFound = table2.reverseLookupX(table2_ls.lookup(xTest, yTest), yTest);
    yFound = table2_ls.reverseLookupY(table2_ls.lookup(xTest, yTest), xTest);
    print("f", table2_ls.lookup(xTest, yTest), func2(xTest, yTest));
    print("dfdx", table2_ls.dFdX(xTest, yTest), dfunc2dx(xTest, yTest));
    print("dfdy", table2_ls.dFdY(xTest, yTest), dfunc2dy(xTest, yTest));
    print("d2fdx2", table2_ls.d2FdX2(xTest, yTest), d2func2dx2(xTest, yTest));
    print("d2fdy2", table2_ls.d2FdY2(xTest, yTest), d2func2dy2(xTest, yTest));
    print("d2fdxdy", table2_ls.d2FdXdY(xTest, yTest), d2func2dxdy(xTest, yTest));
    print("reverseX", xFound, xTest);
    print("reverseY", yFound, yTest);

    Info<<nl<<"3D table" << endl;
    scalarLookupTable3D table3(dict.subDict("table3D"), "x", "y", "z", "f");
    xFound = table3.reverseLookupX(table3.lookup(xTest, yTest, zTest), yTest, zTest);
    yFound = table3.reverseLookupY(table3.lookup(xTest, yTest, zTest), xTest, zTest);
    scalar zFound = table3.reverseLookupZ(table3.lookup(xTest, yTest, zTest), xTest, yTest);
    print("f", table3.lookup(xTest, yTest, zTest), func3(xTest, yTest, zTest));
    print("dfdx", table3.dFdX(xTest, yTest, zTest), dfunc3dx(xTest, yTest, zTest));
    print("dfdy", table3.dFdY(xTest, yTest, zTest), dfunc3dy(xTest, yTest, zTest));
    print("dfdz", table3.dFdZ(xTest, yTest, zTest), dfunc3dz(xTest, yTest, zTest));
    print("d2fdx2", table3.d2FdX2(xTest, yTest, zTest), d2func3dx2(xTest, yTest, zTest));
    print("d2fdy2", table3.d2FdY2(xTest, yTest, zTest), d2func3dy2(xTest, yTest, zTest));
    print("d2fdz2", table3.d2FdZ2(xTest, yTest, zTest), d2func3dz2(xTest, yTest, zTest));
    print("d2fdxdy", table3.d2FdXdY(xTest, yTest, zTest), d2func3dxdy(xTest, yTest, zTest));
    print("d2fdxdz", table3.d2FdXdZ(xTest, yTest, zTest), d2func3dxdz(xTest, yTest, zTest));
    print("d2fdydz", table3.d2FdYdZ(xTest, yTest, zTest), d2func3dydz(xTest, yTest, zTest));
    print("reverseX", xFound, xTest);
    print("reverseY", yFound, yTest);
    print("reverseZ", zFound, zTest);

    Info<<nl<<"3D table from least squares" << endl;
    scalarLookupTable3D table3_ls(dict.subDict("table3D_ls"), "x", "y", "z", "f");
    xFound = table3_ls.reverseLookupX(table3_ls.lookup(xTest, yTest, zTest), yTest, zTest);
    yFound = table3_ls.reverseLookupY(table3_ls.lookup(xTest, yTest, zTest), xTest, zTest);
    zFound = table3_ls.reverseLookupZ(table3_ls.lookup(xTest, yTest, zTest), xTest, yTest);
    print("f", table3_ls.lookup(xTest, yTest, zTest), func3(xTest, yTest, zTest));
    print("dfdx", table3_ls.dFdX(xTest, yTest, zTest), dfunc3dx(xTest, yTest, zTest));
    print("dfdy", table3_ls.dFdY(xTest, yTest, zTest), dfunc3dy(xTest, yTest, zTest));
    print("dfdz", table3_ls.dFdZ(xTest, yTest, zTest), dfunc3dz(xTest, yTest, zTest));
    print("d2fdx2", table3_ls.d2FdX2(xTest, yTest, zTest), d2func3dx2(xTest, yTest, zTest));
    print("d2fdy2", table3_ls.d2FdY2(xTest, yTest, zTest), d2func3dy2(xTest, yTest, zTest));
    print("d2fdz2", table3_ls.d2FdZ2(xTest, yTest, zTest), d2func3dz2(xTest, yTest, zTest));
    print("d2fdxdy", table3_ls.d2FdXdY(xTest, yTest, zTest), d2func3dxdy(xTest, yTest, zTest));
    print("d2fdxdz", table3_ls.d2FdXdZ(xTest, yTest, zTest), d2func3dxdz(xTest, yTest, zTest));
    print("d2fdydz", table3_ls.d2FdYdZ(xTest, yTest, zTest), d2func3dydz(xTest, yTest, zTest));
    print("reverseX", xFound, xTest);
    print("reverseY", yFound, yTest);
    print("reverseZ", zFound, zTest);

    Info<< nl << "Finished" << nl << endl;
    return 0;
}
