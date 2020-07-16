#include "readGravitationalAcceleration.H"

Info<< "Creating twoPhaseSystem\n" << endl;
phaseSystem fluid(mesh);
autoPtr<timeIntegrator> integrator(timeIntegrator::New(mesh));
integrator->addSystem(fluid);

PtrList<phaseModel>& phases(fluid.phases());
