#include "readGravitationalAcceleration.H"

Info<< "Creating twoPhaseSystem\n" << endl;
autoPtr<timeIntegrator> integrator(timeIntegrator::New(mesh));
phaseSystem fluid(mesh);
integrator->addSystem(fluid);

PtrList<phaseModel>& phases(fluid.phases());
