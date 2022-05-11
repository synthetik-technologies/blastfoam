#include "fvCFD.H"
#include "fvMesh.H"
#include "extendedNLevelGlobalCellToCellStencils.H"
#include "Random.H"

using namespace Foam;

template<class StencilType>
void setCells
(
    const label gCelli,
    const label nLevels,
    const StencilType& stencil,
    volScalarField& levels
)
{
    levels == -1.0;
    const Map<cellStencil>& cellCells(stencil.cellCellMap());
    if (cellCells.found(gCelli))
    {
    	const cellStencil& cs = cellCells[gCelli];
    	for (label leveli = 0; leveli <= nLevels; leveli++)
    	{
            const SubList<label> ls = cs.localStencil(leveli);
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
	    do
        {
	       gCelli = ranGen.sampleAB(label(0), nCells);
	    } while (used.found(gCelli));
    	reduce(gCelli, maxOp<label>());
    	samples[i] = gCelli;
    	used.insert(gCelli);
    }

    scalar t0 = runTime.elapsedCpuTime();
    scalar t = t0;
    for (label nLevels = 0; nLevels < 6; nLevels++)
    {
    	Info<< nl << nl << nLevels << " levels" << nl << endl;
    	runTime.setTime(1, 1);
        	extendedNLevelCFCCellToCellStencil stencilF(mesh, nLevels);
    	forAll(samples, i)
    	{
    	    setCells(samples[i], nLevels, stencilF, CFC);
    	    CFC.write();
    	    runTime++;
    	}
    	t = runTime.elapsedCpuTime();
    	Info<< "Face stencil: " << t - t0 << " s" << endl;
    	t0 = t;


    	runTime.setTime(1, 1);
    	extendedNLevelCECCellToCellStencil stencilE(mesh, nLevels);
    	forAll(samples, i)
    	{
    	    setCells(samples[i], nLevels, stencilE, CEC);
    	    CEC.write();
    	    runTime++;
    	}

    	t = runTime.elapsedCpuTime();
    	Info<< "Edge stencil: " << t - t0 << " s" << endl;
    	t0 = t;

    	runTime.setTime(1, 1);
    	extendedNLevelCPCCellToCellStencil stencilP(mesh, nLevels);
    	forAll(samples, i)
    	{
    	    setCells(samples[i], nLevels, stencilP, CPC);
    	    CPC.write();
    	    runTime++;
    	}
    	t = runTime.elapsedCpuTime();
    	Info<< "Point stencil: " << t - t0 << " s" << endl;
    	t0 = t;
    }

    Info<< nl << "done" <<endl;
    return 0;
}
