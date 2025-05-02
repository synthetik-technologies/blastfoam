#include "readGravitationalAcceleration.H"

Info<< "Creating phaseSystem\n" << endl;
fvTimeIntegrator integrator(mesh);
phaseSystem fluid(mesh);
integrator.addSystem(fluid);

PtrList<phaseModel>& phases(fluid.phases());

const ::Foam::fvModels& fvModels = integrator.models();
const ::Foam::fvConstraints& fvConstraints = integrator.constraints();
