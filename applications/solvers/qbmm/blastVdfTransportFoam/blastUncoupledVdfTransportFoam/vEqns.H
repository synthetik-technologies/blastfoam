// Solve the change in velocity momentums due to drag
for (label nodei = 0; nodei < nodes.size(); nodei++)
{
    volScalarField alphaRho(nodes[nodei].primaryWeight()*rhop);
    const volVectorField& Ui = nodes[nodei].velocityAbscissae();
    volScalarField d(nodes[nodei].primaryAbscissae()[sizeIndex]);
    d.max(1e-10);
    volScalarField Re
    (
        "Re",
        mag(U - nodes[nodei].velocityAbscissae())*d/fluid->thermo().nu()
    );
    tmp<volScalarField> Kd;

    if (dragModel == "SchillerNaumann")
    {
        tmp<volScalarField> CdRe
        (
            neg(Re - 1000)*24.0*(1.0 + 0.15*pow(Re, 0.687))
          + pos0(Re - 1000)*0.44*max(Re, 1e-3)
        );
        Kd =
            max(nodes[nodei].primaryWeight(), residualAlpha)
           *0.75*CdRe*rho*fluid->thermo().nu()/sqr(d);

    }
    else if (dragModel == "Stokes")
    {
        dimensionedScalar dragTime
        (
            "dragTime",
            dimTime,
            populationBalanceProperties.subDict("dragCoeffs")
        );
        Kd = max(nodes[nodei].primaryWeight(), residualAlpha)*rhop/dragTime;
    }
    else
    {
         FatalError
            << args.executable() << " : Unknown drag model " << dragModel
            << abort(FatalError);
    }

    fvVectorMatrix UEqn
    (
        alphaRho*fvm::ddt(Ui)
      - alphaRho*fvc::ddt(Ui)
     ==
        Kd()*U
      - fvm::Sp(Kd(), Ui)
      + g*alphaRho
    );
    UEqn.relax();
    UEqn.solve();
}
quadrature.updateMoments();
