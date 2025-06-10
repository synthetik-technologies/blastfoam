fvTimeIntegrator integrator(mesh);
autoPtr<compressibleSystem> fluid
(
    compressibleSystem::New(mesh)
);
integrator.addSystem(fluid());
fluid->update();

const volScalarField& p = fluid->p();
const volScalarField& T = fluid->T();

const fvModels& fvModels = integrator.models();
const fvConstraints& fvConstraints = integrator.constraints();
