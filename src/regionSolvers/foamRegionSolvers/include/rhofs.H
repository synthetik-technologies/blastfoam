#ifdef COMPRESSIBLE
surfaceScalarField rho1f(fvc::interpolate(rho1));
surfaceScalarField rho2f(fvc::interpolate(rho2));
#else
const dimensionedScalar& rho1f(rho1);
const dimensionedScalar& rho2f(rho2);
#endif