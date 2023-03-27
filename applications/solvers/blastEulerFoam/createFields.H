#include "readGravitationalAcceleration.H"

Info<< "Creating phaseSystem\n" << endl;
fvTimeIntegrator integrator(mesh);
phaseSystem fluid(mesh);
integrator.addSystem(fluid);

PtrList<phaseModel>& phases(fluid.phases());
