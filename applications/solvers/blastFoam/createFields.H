fvTimeIntegrator integrator(mesh);
autoPtr<compressibleSystem> fluid
(
    compressibleSystem::New(mesh)
);
integrator.addSystem(fluid());

const volScalarField& p = fluid->p();
const volScalarField& T = fluid->T();

integrator.models();
integrator.constraints();

// fluid->update();
// integrator.clear();
