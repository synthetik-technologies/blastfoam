// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> deltaUEqn
(
    fvm::div(phi, U)
  + turbulence->divDevTau(U)
);
tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U)
  + MRF.DDt(rho, U)
  + deltaUEqn()
 ==
    fvModels.source(rho, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();
fvConstraints.constrain(UEqn);

volVectorField gradp
(
    -fvc::reconstruct((fvc::snGrad(p))*mesh.magSf())
);

// solve(UEqn == gradp);
// phi = fvc::flux(U*rho);

volVectorField F
(
    IBM.forcing
    (
        U,
        rho,
        (rho.oldTime()*U.oldTime())(),
        ((deltaUEqn & U) + fvc::grad(p))(),
        runTime.deltaT()
    )
);

UEqn -= F;

UEqn.relax();
fvConstraints.constrain(U);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        gradp
    );

    fvConstraints.constrain(U);
    K = 0.5*magSqr(U);
}
