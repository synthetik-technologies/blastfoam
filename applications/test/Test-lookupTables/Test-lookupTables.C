#include "dictionary.H"
#include "scalarLookupTable1D.H"
#include "scalarLookupTable2D.H"
#include "scalarLookupTable3D.H"
#include "OFstream.H"
#include "IFstream.H"
#include "argList.H"

using namespace Foam;

scalar func(const scalar x, const scalar y, const scalar z)
{
    return 3.0*x + y + sqr(z);
}

int main(int argc, char *argv[])
{
    // Create some tables
    label nx = 5;
    label ny = 7;
    label nz = 9;
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

//     {
//         OFstream out1D("table1D.csv");
//         for (label i = 0; i < nx; i++)
//         {
//            out1D << 2*(x[i]);
//             if ( i != nx-1)
//             {
//                 out1D << ",";
//             }
//         }
//     }
//     {
//         OFstream out2D("table2D.csv");
//         for (label j = 0; j < ny; j++)
//         {
//             for (label i = 0; i < nx; i++)
//             {
//                 out2D << (x[i]) + (y[j]);
//                 if ( i != nx-1)
//                 {
//                     out2D << ",";
//                 }
//             }
//             out2D <<endl;
//         }
//     }
    {
        OFstream out3D("table3D.csv");
        for (label k = 0; k < nz; k++)
        {
            for (label j = 0; j < ny; j++)
            {
                for (label i = 0; i < nx; i++)
                {
                    out3D << func(x[i], y[j], z[k]);
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

    scalar p = 1e6;
    scalar T = 1500;

    Info<<nl<<"1D table:" << endl;
    scalarLookupTable1D table1(dict.subDict("table1D"), "T", "e");
    scalar e = table1.lookup(T);
    Info<< "e: " << e <<endl;
    Info<< "T: " << table1.reverseLookup(e) <<endl;
    Info<< nl << nl;

    Info<<nl<<"2D table:" << endl;
    scalarLookupTable2D table2(dict.subDict("table2D"), "rho", "e", "p");
    scalar rho = table2.reverseLookupX(p, e);
    Info<< "rho: " << rho <<endl;
    Info<< "p: " << table2.lookup(rho, e) <<endl;
    Info<< "T: " << table1.reverseLookup(table2.reverseLookupY(p, rho)) <<endl;

    Info<<nl<<"3D table" << endl;
    scalarLookupTable3D table3(dict.subDict("table3D"), "x", "y", "z", "f");
    scalar xTest = 0.1;
    scalar yTest = 1.6;
    scalar zTest = 2.6;

    scalar f = table3.lookup(xTest, yTest, zTest);
    Info<< "f: " << f <<endl;
    Info<< "answer: " << func(xTest, yTest, zTest) << endl;

    return 0;
}
