#include "fvCFD.H"
#include "dynamicBlastFvMesh.H"
#include "immersedBoundaryObjectListSolver.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicBlastFvMesh.H"

    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedVector(dimAcceleration, Zero)
    );
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPressure, 0)
    );

    immersedBoundaryObjectListSolver ibm(mesh);

    while (runTime.run())
    {
        mesh.refine();
        Info<<"Time: " << runTime.timeName() << endl;
        ibm.solve();
        runTime.write();
        runTime++;
    }
    runTime.write();
    Info<<"done"<<endl;

    return 0;
}
