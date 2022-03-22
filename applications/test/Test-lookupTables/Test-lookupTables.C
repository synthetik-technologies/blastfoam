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

scalar func2(const scalar x, const scalar y)
{
    return 3.0*x + sqr(y)*x;
}

scalar func3(const scalar x, const scalar y, const scalar z)
{
    return 3.0*x*Foam::sin(y) + y*Foam::pow(z, 1.3);
}

int main(int argc, char *argv[])
{
    // Create some tables
    label nx = 15;
    label ny = 27;
    label nz = 39;
    scalar xMax = 2.0;
    scalar yMax = 3.0;
    scalar zMax = 4.0;
    scalar dx = xMax/scalar(nx-1.0);
    scalar dy = yMax/scalar(ny-1.0);
    scalar dz = zMax/scalar(nz-1.0);

    scalarField x(nx);
    scalarField y(ny);
    scalarField z(nz);

    {
        OFstream outX("x.csv");
        for (label i = 0; i < nx; i++)
        {
            x[i] =  dx*scalar(i);
            outX << x[i] << ";";
        }
    }
    {
        OFstream outY("y.csv");
        for (label j = 0; j < ny; j++)
        {
            y[j] =  dy*scalar(j);
            outY << y[j] << ";\n";
        }
    }
    {
        OFstream outZ("z.csv");
        outZ << "# abc"<<endl;
        for (label k = 0; k < nz; k++)
        {
            z[k] =  dz*scalar(k);
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

    scalar xTest = 1.5;
    scalar yTest = 1.6;
    scalar zTest = 2.6;

    Info<<nl<<"1D table:" << endl;
    scalarLookupTable1D table1(dict.subDict("table1D"), "x", "f");
    scalarLookupTable1D table11(table1);
    Info<< "f: " << table1.lookup(xTest) <<endl;
    Info<< "fCopy: " << table11.lookup(xTest) <<endl;
    Info<< "answer: " << func1(xTest) << endl;

    Info<<nl<<"2D table:" << endl;
    lookupTable2D<scalar> table2(dict.subDict("table2D"), "x", "y", "f");
    scalarLookupTable2D table21(table2);
    Info<< "f: " << table2.lookup(xTest, yTest) <<endl;
    Info<< "f interpolate: " << table2.interpolate(table2.f()) <<endl;
    Info<< "fCopy: " << table21.lookup(xTest, yTest) <<endl;
    Info<< "answer: " << func2(xTest, yTest) << endl;

    Info<<nl<<"3D table" << endl;
    scalarLookupTable3D table3(dict.subDict("table3D"), "x", "y", "z", "f");
    scalarLookupTable3D table31(table3);
    Info<< "f: " << table3.lookup(xTest, yTest, zTest) <<endl;
    Info<< "f interpolate: " << table3.interpolate(table3.f()) <<endl;
    Info<< "fCopy: " << table31.lookup(xTest, yTest, zTest) <<endl;
    Info<< "answer: " << func3(xTest, yTest, zTest) << endl;

    return 0;
}
