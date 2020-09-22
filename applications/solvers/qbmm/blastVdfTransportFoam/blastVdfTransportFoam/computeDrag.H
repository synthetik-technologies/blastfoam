// Solve the change in velocity momentums due to drag
PtrList<volScalarField> Kds(nodes.size());
for (label nodei = 0; nodei < nodes.size(); nodei++)
{
    const volVectorField& Ui = nodes[nodei].velocityAbscissae();
    volScalarField d(nodes[nodei].primaryAbscissae()[sizeIndex]);
    d.max(1e-8);
    volScalarField Re
    (
        "Re",
        max(mag(U - Ui)*d/fluid.thermo().nu(), 1e-3)
    );
    const volScalarField& alphai = nodes[nodei].primaryWeight();

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

        if (dragModel == "WenYu" || dragModel == "GidaspowErgunWenYu")
        {
            volScalarField alpha2(max(alphac, residualAlpha));
            volScalarField Res(alpha2*Re);
            volScalarField CdsRes
            (
                neg(Res - 1000)*24.0*(1.0 + 0.15*pow(Res, 0.687))
              + pos0(Res - 1000)*0.44*max(Res, 1e-3)
            );

            CdReWenYu =
            (
                CdsRes
               *pow(alpha2, -3.65)
               *max(alphac, residualAlpha)
            );
        }
        if (dragModel == "Ergun" || dragModel == "GidaspowErgunWenYu")
        {
            CdReErgun =
                (4.0/3.0)
               *(
                    150*(1.0 - alphac)/max(alphac, residualAlpha)
                  + 1.75*Re
                );
        }

        if (CdReWenYu.valid() && CdReErgun.valid())
        {
            CdRe =
                pos0(alphac - 0.8)*CdReWenYu
              + neg(alphac - 0.8)*CdReErgun;
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

    volScalarField Kdi(0.75*CdRe*fluid.thermo().mu()/sqr(d));
    Kds.set
    (
        nodei,
        max(alphai, residualAlpha)*Kdi
    );

    fluid.addUCoeff(Kds[nodei]);
    fluid.addUSource(Kds[nodei]*Ui);
}
