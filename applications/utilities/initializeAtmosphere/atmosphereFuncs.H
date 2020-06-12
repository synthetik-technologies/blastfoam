namespace Foam
{

void tableAtmosphere
(
    const dictionary& dict,
    const uniformDimensionedVectorField& g,
    const vector& dir,
    const volScalarField& h,
    volScalarField& p,
    volScalarField& rho
)
{
    dictionary pTableDict(dict.subDict("pTable"));
    lookupTable1D pTable
    (
        pTableDict.lookup("file"),
        pTableDict.lookupOrDefault<word>("pMod", "none"),
        pTableDict.lookupOrDefault<word>("hMod", "none")
    );
    dictionary TTableDict(dict.subDict("TTable"));
    lookupTable1D TTable
    (
        pTableDict.lookup("file"),
        pTableDict.lookupOrDefault<word>("TMod", "none"),
        pTableDict.lookupOrDefault<word>("hMod", "none")
    );

    scalar R = 8314.32/28.97;
    forAll(rho, celli)
    {
        p[celli] = pTable.lookup(h[celli]);
        scalar T = TTable.lookup(h[celli]);
        rho[celli] = p[celli]/T/R;
    }

    forAll(rho.boundaryField(), patchi)
    {
        forAll(rho.boundaryField()[patchi], facei)
        {
            p.boundaryFieldRef()[patchi][facei] =
                pTable.lookup(h.boundaryField()[patchi][facei]);
            scalar T = TTable.lookup(h.boundaryField()[patchi][facei]);
            rho.boundaryFieldRef()[patchi][facei] =
                p.boundaryField()[patchi][facei]/T/R;
        }
    }
}

void simpleAtmosphere
(
    const dictionary& dict,
    const uniformDimensionedVectorField& g,
    const vector& dir,
    const volScalarField& h,
    volScalarField& p,
    volScalarField& rho
)
{
    dimensionedScalar pRef("pRef", dimPressure, dict);
    p = pRef + (g & dir)*rho*h;
    p.correctBoundaryConditions();
}

}
