#include "dictionary.H"
#include "lookupTables1D.H"
#include "lookupTables2D.H"
#include "lookupTables3D.H"
#include "OFstream.H"
#include "IFstream.H"
#include "argList.H"

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

int main(int argc, char *argv[])
{
    // Create some tables
    label nx = 20;
    label ny = 20;
    label nz = 30;
    scalar xMin = 1.0;
    scalar yMin = 0.1;
    scalar zMin = 0.001;
    scalar xMax = 2.0;
    scalar yMax = 3.0;
    scalar zMax = 4.0;
    scalar dx = (xMax - xMin)/scalar(nx);
    scalar dy = (yMax - yMin)/scalar(ny);
    scalar dz = (zMax - zMin)/scalar(nz);

    scalarField x(nx);
    scalarField y(ny);
    scalarField z(nz);

    {
        OFstream outX("x.csv");
        for (label i = 0; i < nx; i++)
        {
            x[i] =  xMin + dx*scalar(i);
            outX << x[i] << ";";
        }
    }
    {
        OFstream outY("y.csv");
        for (label j = 0; j < ny; j++)
        {
            y[j] =  yMin + dy*scalar(j);
            outY << y[j] << ";\n";
        }
    }
    {
        OFstream outZ("z.csv");
        outZ << "# abc"<<endl;
        for (label k = 0; k < nz; k++)
        {
            z[k] =  zMin + dz*scalar(k);
            outZ << z[k] << "\n";
        }
    }

    {
        OFstream out1D("table1D.csv");
        for (label i = 0; i < nx; i++)
        {
            out1D << func1(x[i]);
            if ( i != nx-1)
            {
                out1D << ",";
            }
        }
    }
    {
        OFstream out2D("table2D.csv");
        for (label j = 0; j < ny; j++)
        {
            for (label i = 0; i < nx; i++)
            {
                out2D << func2(x[i], y[j]);
                if ( i != nx-1)
                {
                    out2D << ",";
                }
            }
            out2D <<endl;
        }
    }
    {
        OFstream out3D("table3D.csv");
        for (label k = 0; k < nz; k++)
        {
            for (label j = 0; j < ny; j++)
            {
                for (label i = 0; i < nx; i++)
                {
                    out3D << func3(x[i], y[j], z[k]);
                    if ( i != nx-1)
                    {
                        out3D << ",";
                    }
                }
                out3D <<";";
            }
            out3D<<endl;
        }
    }


    IFstream is("tableDict");
    dictionary dict(is);

    scalar xTest = 1.435;
    scalar yTest = 1.3346;
    scalar zTest = 2.5676;

    Info<<nl<<"1D table:" << endl;
    scalarLookupTable1D table1(dict.subDict("table1D"), "x", "f");
    scalarLookupTable1D table11(table1);
    Info<< "f: " << table1.lookup(xTest)
        << ", answer: " << func1(xTest) << endl
        << "dfdx: " << table1.dFdX(xTest)
        << ", answer: " << dfunc1dx(xTest) << endl
        << "d2fdx2: " << table1.d2FdX2(xTest)
        << ", answer: " << d2func1dx2(xTest) << endl;

    Info<<nl<<"2D table:" << endl;
    lookupTable2D<scalar> table2(dict.subDict("table2D"), "x", "y", "f");
    scalarLookupTable2D table21(table2);
    Info<< "f: " << table2.lookup(xTest, yTest)
        << ", answer: " << func2(xTest, yTest) << endl
        << "dfdx: " << table2.dFdX(xTest, yTest)
        << ", answer: " << dfunc2dx(xTest, yTest) << endl
        << "dfdy: " << table2.dFdY(xTest, yTest)
        << ", answer: " << dfunc2dy(xTest, yTest) << endl
        << "d2fdx2: " << table2.d2FdX2(xTest, yTest)
        << ", answer: " << d2func2dx2(xTest, yTest) << endl
        << "d2fdy2: " << table2.d2FdY2(xTest, yTest)
        << ", answer: " << d2func2dy2(xTest, yTest) << endl
        << "d2fdxdy: " << table2.d2FdXdY(xTest, yTest)
        << ", answer: " << d2func2dxdy(xTest, yTest) << endl
        << "reverseX: "<< table2.reverseLookupX(table2.lookup(xTest, yTest), yTest)<< ", answer: "<<xTest<<endl
        << "reverseY: "<< table2.reverseLookupY(table2.lookup(xTest, yTest), xTest)<< ", answer: "<<yTest<<endl;

    Info<<nl<<"3D table" << endl;
    scalarLookupTable3D table3(dict.subDict("table3D"), "x", "y", "z", "f");
    scalarLookupTable3D table31(table3);
    Info<< "f: " << table3.lookup(xTest, yTest, zTest)
        << ", answer: " << func3(xTest, yTest, zTest) << endl
        << "dfdx: " << table3.dFdX(xTest, yTest, zTest)
        << ", answer: " << dfunc3dx(xTest, yTest, zTest) << endl
        << "dfdy: " << table3.dFdY(xTest, yTest, zTest)
        << ", answer: " << dfunc3dy(xTest, yTest, zTest) << endl
        << "dfdz: " << table3.dFdZ(xTest, yTest, zTest)
        << ", answer: " << dfunc3dz(xTest, yTest, zTest) << endl
        << "d2fdx2: " << table3.d2FdX2(xTest, yTest, zTest)
        << ", answer: " << d2func3dx2(xTest, yTest, zTest) << endl
        << "d2fdy2: " << table3.d2FdY2(xTest, yTest, zTest)
        << ", answer: " << d2func3dy2(xTest, yTest, zTest) << endl
        << "d2fdz2: " << table3.d2FdZ2(xTest, yTest, zTest)
        << ", answer: " << d2func3dz2(xTest, yTest, zTest) << endl
        << "d2fdxdy: " << table3.d2FdXdY(xTest, yTest, zTest)
        << ", answer: " << d2func3dxdy(xTest, yTest, zTest) << endl
        << "d2fdxdz: " << table3.d2FdXdZ(xTest, yTest, zTest)
        << ", answer: " << d2func3dxdz(xTest, yTest, zTest) << endl
        << "d2fdydz: " << table3.d2FdYdZ(xTest, yTest, zTest)
        << ", answer: " << d2func3dydz(xTest, yTest, zTest) << endl;

    Info<< nl << "Finished" << nl << endl;
    return 0;
}
