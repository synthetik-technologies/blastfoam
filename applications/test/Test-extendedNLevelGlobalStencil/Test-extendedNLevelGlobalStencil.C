#include "fvCFD.H"
#include "fvMesh.H"
#include "extendedNLevelGlobalCellToCellStencils.H"
#include "Random.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    runTime.setDeltaT(1);
    volScalarField CFC
    (
        IOobject
        (
            "CFC",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        -1.0
    );
    volScalarField CEC
    (
        IOobject
        (
            "CEC",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        -1.0
    );
    volScalarField CPC
    (
        IOobject
        (
            "CPC",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        -1.0
    );

    Random ranGen(123);
    const label nCells = returnReduce(mesh.nCells(), sumOp<label>());
    const label nSamples = min(25, nCells);

    const label nLevels = 3;
    extendedNLevelCFCCellToCellStencil stencilF(mesh, nLevels);
    extendedNLevelCECCellToCellStencil stencilE(mesh, nLevels);
    extendedNLevelCPCCellToCellStencil stencilP(mesh, nLevels);

    const globalIndex gIndex(mesh.nCells());
    labelHashSet used;
    label gCelli = 0;
    do
    {
        if (Pstream::master())
        {
            while (used.found(gCelli))
            {
                gCelli = ranGen.sampleAB(label(0), nCells);
            }
        }
        else
        {
            gCelli = -1;
        }
        reduce(gCelli, maxOp<label>());
        used.insert(gCelli);

        CFC == -1.0;
        CEC == -1.0;
        CPC == -1.0;

        for (label leveli = 0; leveli <= nLevels; leveli++)
        {
            const Map<cellStencil>& cellCells(stencilF.cellCellMap());
            if (cellCells.found(gCelli))
            {
                const cellStencil& cs = cellCells[gCelli];
                const labelList& ls = cs.localStencil(leveli);
                forAll(ls, i)
                {
                    CFC[ls[i]] = leveli;
                }
            }
        }
        for (label leveli = 0; leveli <= nLevels; leveli++)
        {
            const Map<cellStencil>& cellCells(stencilE.cellCellMap());
            if (cellCells.found(gCelli))
            {
                const cellStencil& cs = cellCells[gCelli];
                const labelList& ls = cs.localStencil(leveli);
                forAll(ls, i)
                {
                    CEC[ls[i]] = leveli;
                }
            }
        }
        for (label leveli = 0; leveli <= nLevels; leveli++)
        {
            const Map<cellStencil>& cellCells(stencilP.cellCellMap());
            if (cellCells.found(gCelli))
            {
                const cellStencil& cs = cellCells[gCelli];
                const labelList& ls = cs.localStencil(leveli);
                forAll(ls, i)
                {
                    CPC[ls[i]] = leveli;
                }
            }
        }

        runTime.write();
        runTime++;
    } while (runTime.timeIndex() < nSamples);

    Info<<"done"<<endl;

    return 0;
}
