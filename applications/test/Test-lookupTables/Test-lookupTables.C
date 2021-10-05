#include "dictionary.H"
#include "scalarLookupTable1D.H"
#include "scalarLookupTable2D.H"
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
    return 0;
}
