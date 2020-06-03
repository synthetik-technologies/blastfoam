// if (!isA<staticFvMesh>(mesh))
// {
//     error.update();
//     mesh.update();
// }

fluid.encode();

Info<< "Calculating Fluxes" << endl;
integrator.integrate();

Info<< "max(p): " << max(fluid.p()).value()
    << ", min(p): " << min(fluid.p()).value() << endl;
Info<< "max(T): " << max(fluid.T()).value()
    << ", min(T): " << min(fluid.T()).value() << endl;

fluid.clearODEFields();
