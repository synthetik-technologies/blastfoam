#include "dictionary.H"
#include "scalarLookupTable1D.H"
#include "scalarLookupTable2D.H"
#include "scalarLookupTable3D.H"
#include "OFstream.H"
#include "IFstream.H"
#include "argList.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    IFstream is("tableDict");
    dictionary dict(is);

    scalar p = 1e6;
    scalar T = 1500;

    Info<<"1D table:" << endl;
    scalarLookupTable1D table1(dict.subDict("table1D"), "T", "e");
    scalar e = 1.63468e+06;//table1.lookup(T);
    Info<< "e: " << e <<endl;
    Info<< "T: " << table1.reverseLookup(e) <<endl;
    Info<< nl << nl;

    Info<<"2D table:" << endl;
    scalarLookupTable2D table2(dict.subDict("table2D"), "rho", "e", "p");
    scalar rho = table2.reverseLookupX(p, e);
    Info<< "rho: " << rho <<endl;
    Info<< "p: " << table2.lookup(rho, e) <<endl;
    Info<< "T: " << table1.reverseLookup(table2.reverseLookupY(p, rho)) <<endl;

//     OFstream out3D("table3D.csv");
//     label nx = 5;
//     label ny = 7;
//     label nz = 9;
//     scalar xMax = 2.0;
//     scalar yMax = 3.0;
//     scalar zMax = 4.0;
//     scalar dx = xMax/scalar(nx-1.0);
//     scalar dy = yMax/scalar(ny-1.0);
//     scalar dz = zMax/scalar(nz-1.0);
//     for (label k = 0; k < nz; k++)
//     {
//         scalar z =  dz*scalar(k);
//         for (label j = 0; j < ny; j++)
//         {
//             scalar y =  dy*scalar(j);
//             for (label i = 0; i < nx; i++)
//             {
//                 scalar x =  dx*scalar(i);
//                 out3D << 3*(x) + (y) + (z);
//                 if ( i != nx-1)
//                 {
//                     out3D << ",";
//                 }
//             }
//             out3D <<";";
//         }
//         out3D<<endl;
//     }
    Info<<"3D table:" << endl;
    scalarLookupTable3D table3(dict.subDict("table3D"), "x", "y", "z", "f");
    scalar f = table3.lookup(0.1, 1.6, 2.6);
    Info<< "f: " << f <<endl;

    return 0;
}
