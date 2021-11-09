#include "fvCFD.H"
#include "fvMesh.H"
#include "immersedBoundaryObjectListSolver.H"

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
    volScalarField internalExternal
    (
        IOobject
        (
            "internalExternal",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        0.0
    );
    volScalarField shell
    (
        IOobject
        (
            "shell",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        0.0
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
    immersedBoundaryObjectListSolver& ibm =
        immersedBoundaryObjectListSolver::New(mesh);
    const PtrList<immersedBoundaryObject>& objects(ibm.objects());
    ibm.setCellTypes();
    ibm.setObjectIDs();

    forAll(objects, i)
    {
        const labelList& pI = objects[i].patchInternalCells();
        const labelList& pE = objects[i].patchExternalCells();
        forAll(pI, j)
        {
            label celli = pI[j];
            if (celli > -1)
            {
                internalExternal[celli] = 1.0;
            }
        }
        forAll(pE, j)
        {
            label celli = pE[j];
            if (celli > -1)
            {
                internalExternal[celli] = 2.0;
            }
        }
        objects[i].setShell(shell, 1.0);
        objects[i].writeVTK(objects[i].name());
    }

//     runTime++;
    runTime.writeNow();
    Info<<"done"<<endl;

    return 0;
}
