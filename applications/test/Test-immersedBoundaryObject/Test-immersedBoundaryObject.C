#include "fvCFD.H"
#include "fvMesh.H"
#include "immersedBoundaryObjects.H"
#include "emptyPolyPatch.H"
#include "wedgePolyPatch.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    volScalarField cellType
    (
        IOobject
        (
            "cellType",
            runTime.timeName(),
            mesh
        ),
        mesh,
        -1
    );
    volVectorField normals
    (
        IOobject
        (
            "normals",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedVector(dimless, Zero)
    );

    IOdictionary ibmDict
    (
        IOobject
        (
            "immersedBoundaryProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ
        )
    );

    const dictionary& objectDict(ibmDict.subDict("objects"));
    wordList objects(objectDict.toc());

    PtrList<immersedBoundaryObject> ibmObjects(objects.size());
    forAll(objects, i)
    {
        ibmObjects.set
        (
            i,
            immersedBoundaryObject::New
            (
                mesh,
                objectDict.subDict(objects[i]),
                objectDict.subDict(objects[i])
            ).ptr()
        );
        ibmObjects[i].initialize();
        ibmObjects[i].setInternal(cellType, 1.0, maxEqOp<scalar>());
        ibmObjects[i].setShell(cellType, 2.0, maxEqOp<scalar>());
        ibmObjects[i].setBoundary(cellType, 3.0, maxEqOp<scalar>());
//         labelList pI(ibmObjects[i].patchInternalCells());
//         labelList pE(ibmObjects[i].patchExternalCells());
//         forAll(pI, j)
//         {
//             if (pI[j] >= 0)
//             {
//                 cellType[pI[j]] = 1.0;
//             }
//         }
//         forAll(pE, j)
//         {
//             if (pE[j] >= 0)
//             {
//                 cellType[pE[j]] = 2.0;
//             }
//         }

        ibmObjects[i].shape().writeVTK();

        vectorField n(ibmObjects[i].Sf());
        ibmObjects[i].interpolateFrom(n, normals.primitiveFieldRef());
    }

//     runTime++;
    runTime.write();
    normals.write();
    cellType.write();
    Info<<"done"<<endl;

    return 0;
}
