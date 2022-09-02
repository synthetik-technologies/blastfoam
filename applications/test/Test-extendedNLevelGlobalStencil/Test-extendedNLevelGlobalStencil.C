#include "fvCFD.H"
#include "fvMesh.H"
#include "extendedNLevelGlobalCellToCellStencils.H"
#include "Random.H"

using namespace Foam;

template<class StencilType>
void setCells
(
    const label gCelli,
    const label leveli,
    const StencilType& stencil,
    volScalarField& levels
)
{
    const Map<cellStencil>& cellCells(stencil.cellCellMap());
    if (cellCells.found(gCelli))
    {
	   const cellStencil& cs = cellCells[gCelli];
    	// for (label leveli = nLevels; leveli >= 0; leveli--)
    	{
                const labelList& ls = cs.localStencil();
                forAll(ls, i)
                {
                    levels[ls[i]] = leveli;
                }
        }
    }
}


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

    const globalIndex gIndex(mesh.nCells());

    Random ranGen(123);
    const label nCells = returnReduce(mesh.nCells(), sumOp<label>());
    const label nSamples = min(25, nCells);

    labelList samples(nSamples, -1);
    labelHashSet used;
    forAll(samples, i)
    {
    	label gCelli = -1;
    	if (Pstream::master())
    	{
    	    do
    	    {
                gCelli = ranGen.sampleAB(label(0), nCells);
    	    } while (used.found(gCelli));
    	}
    	reduce(gCelli, maxOp<label>());
    	samples[i] = gCelli;
    	used.insert(gCelli);
    }

    label maxLevels = 5;
    PtrList<extendedNLevelCFCCellToCellStencil> stencilF(maxLevels+1);
    PtrList<extendedNLevelCECCellToCellStencil> stencilE(maxLevels+1);
    PtrList<extendedNLevelCPCCellToCellStencil> stencilP(maxLevels+1);


    for (label nLevels = maxLevels; nLevels > 0; nLevels--)
    {
        Info<< nl << nl << nLevels << " levels" << nl << endl;

        scalar t0 = runTime.elapsedCpuTime();
        scalar t = t0;

        stencilF.set(nLevels, new extendedNLevelCFCCellToCellStencil(mesh, nLevels));
    	t = runTime.elapsedCpuTime();
    	Info<< "Face stencil: " << t - t0 << " s" << endl;
    	t0 = t;


    	stencilE.set(nLevels, new extendedNLevelCECCellToCellStencil(mesh, nLevels));
    	t = runTime.elapsedCpuTime();
    	Info<< "Edge stencil: " << t - t0 << " s" << endl;
    	t0 = t;

    	stencilP.set(nLevels, new extendedNLevelCPCCellToCellStencil(mesh, nLevels));
    	t = runTime.elapsedCpuTime();
    	Info<< "Point stencil: " << t - t0 << " s" << endl;
    	t0 = t;
    }


    forAll(samples, i)
    {
        CFC == -1;
        CEC == -1;
        CPC == -1;
        for (label nLevels = maxLevels; nLevels > 0; nLevels--)
        {
            setCells(samples[i], nLevels, stencilF[nLevels], CFC);
            setCells(samples[i], nLevels, stencilE[nLevels], CEC);
            setCells(samples[i], nLevels, stencilP[nLevels], CPC);
        }

        CFC.write();
        CEC.write();
        CPC.write();
        runTime++;
    }

    Info<< nl << "done" <<endl;
    return 0;
}
