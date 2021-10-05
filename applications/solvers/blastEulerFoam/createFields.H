#include "readGravitationalAcceleration.H"

Info<< "Creating phaseSystem\n" << endl;
autoPtr<timeIntegrator> integrator(timeIntegrator::New(mesh));
phaseSystem fluid(mesh);
integrator->addSystem(fluid);

PtrList<phaseModel>& phases(fluid.phases());
