if (composition.contains("ft"))
{
    volScalarField& ft = composition.Y("ft");

    fvScalarMatrix ftEqn
    (
        fvm::ddt(rho, ft)
      + fvm::div(phi, ft)
      - fvm::laplacian(thermophysicalTransport.alphaEff(), ft)
    );
    ftEqn.solve();
}
