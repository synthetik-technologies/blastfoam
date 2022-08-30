#include "polynomialRoots.H"
#include "OStringStream.H"
#include "GaussianWeights.H"
#include "argList.H"

using namespace Foam;

int main(int argc, char *argv[])
{

    polynomialRoots P1({1.0, 0.0, -1.0});
    List<scalar> rootsReP1
    (
        {1.0, -1.0}
    );
    List<scalar> rootsImP1
    (
        {0.0, 0.0}
    );
    Info<< "f(x) = " << P1.polyName()<<endl;

    Info<< "Roots:" << nl << incrIndent;
    forAll(rootsReP1, i)
    {
        Info<< indent << "real: "<< P1.rootsRe()[i]
            << ", imaginary: " << P1.rootsIm()[i]
            << ", error: "
            << mag(P1.rootsRe()[i] - rootsReP1[i]) + mag(P1.rootsIm()[i] - rootsImP1[i]) << endl;
    }
    Info<< decrIndent << endl;


    polynomialRoots P2({2.34, -2.345, 93.12, 912});
    List<scalar> rootsReP2
    (
        {3.157753198249150, 3.157753198249150, -5.313369644361553}
    );
    List<scalar> rootsImP2
    (
        {7.961161128236527, -7.961161128236527, 0.0}
    );

    Info<< "f(x) = " << P2.polyName()<<endl;
    Info<< "Roots:" << nl << incrIndent;
    forAll(rootsReP2, i)
    {
        Info<< indent << "real: "<< P2.rootsRe()[i]
            << ", imaginary: " << P2.rootsIm()[i]
            << ", error: "
            << mag(P2.rootsRe()[i] - rootsReP2[i]) + mag(P2.rootsIm()[i] - rootsImP2[i]) << endl;
    }
    Info<< decrIndent << endl;


    Info<< nl << "Finished" << nl << endl;
    return 0;
}
