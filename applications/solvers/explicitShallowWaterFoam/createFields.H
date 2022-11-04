autoPtr<timeIntegrator> integrator(timeIntegrator::New(mesh));
shallowWaterSystem fluid(mesh);

integrator->addSystem(fluid);