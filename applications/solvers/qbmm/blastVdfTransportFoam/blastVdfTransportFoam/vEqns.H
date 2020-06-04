// Solve the change in velocity momentums due to drag
for (label nodei = 0; nodei < nodes.size(); nodei++)
{
    volScalarField alphaRho(nodes[nodei].primaryWeight()*rhop);
    const volVectorField& Ui = nodes[nodei].velocityAbscissae();

    fvVectorMatrix UEqn
    (
        alphaRho*fvm::ddt(Ui)
      - alphaRho*fvc::ddt(Ui)
     ==
        Kds[nodei]*U
      - fvm::Sp(Kds[nodei], Ui)
      - nodes[nodei].primaryWeight()*fvc::grad(p)
      + alphaRho*g
    );
    UEqn.relax();
    UEqn.solve();
}
quadrature.updateMoments();
