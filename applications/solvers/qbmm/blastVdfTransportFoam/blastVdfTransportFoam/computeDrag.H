// Solve the change in velocity momentums due to drag
PtrList<volScalarField> Kds(nodes.size());
for (label nodei = 0; nodei < nodes.size(); nodei++)
{
    const volVectorField& Ui = nodes[nodei].velocityAbscissae();
    volScalarField d(nodes[nodei].primaryAbscissae()[sizeIndex]);
    d.max(1e-6);
    volScalarField Re
    (
        "Re",
        mag(U - nodes[nodei].velocityAbscissae())*d/fluid.thermo().nu()
    );
    const volScalarField alphai = nodes[nodei].primaryWeight();

    tmp<volScalarField> CdRe;
    if (dragModel == "SchillerNaumann")
    {
        CdRe =
        (
            neg(Re - 1000)*24.0*(1.0 + 0.15*pow(Re, 0.687))
          + pos0(Re - 1000)*0.44*max(Re, 1e-3)
        );
    }
    else if
    (
        dragModel == "WenYu"
     || dragModel == "Ergun"
     || dragModel == "GidaspowErgunWenYu"
    )
    {
        tmp<volScalarField> CdReWenYu;
        tmp<volScalarField> CdReErgun;
        volScalarField alpha2(max(scalar(1) - alpha, residualAlpha));

        if (dragModel == "WenYu" || dragModel == "GidaspowErgunWenYu")
        {
            volScalarField Res(alpha2*Re);
            volScalarField CdsRes
            (
                neg(Res - 1000)*24.0*(1.0 + 0.15*pow(Res, 0.687))
              + pos0(Res - 1000)*0.44*Res
            );

            CdReWenYu =
            (
                CdsRes
               *pow(alpha2, -2.65)
            );
        }
        if (dragModel == "Ergun" || dragModel == "GidaspowErgunWenYu")
        {
            CdReErgun =
                (4.0/3.0)
               *(
                    150*alphai/max(alpha2, residualAlpha)
                  + 1.75*Re
                );
        }

        if (CdReWenYu.valid() && CdReErgun.valid())
        {
            CdRe =
                pos0(alpha2 - 0.8)*CdReWenYu
              + neg(alpha2 - 0.8)*CdReErgun;
        }
        else if (CdReWenYu.valid())
        {
            CdRe = CdReWenYu;
        }
        else
        {
            CdRe = CdReErgun;
        }
    }
    else
    {
         FatalError
            << args.executable() << " : Unknown drag model " << dragModel
            << abort(FatalError);
    }

    volScalarField Kdi(0.75*CdRe*rho*fluid.thermo().nu()/sqr(d));
    Kds.set
    (
        nodei,
        max(alphai, residualAlpha)*Kdi
    );
    volScalarField KD(Kdi*nodes[nodei].primaryWeight());

    fluid.addUCoeff(KD);
    fluid.addUSource(-KD*Ui);
}
