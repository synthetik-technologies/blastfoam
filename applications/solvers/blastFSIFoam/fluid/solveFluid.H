Info<< "Calculating Fluxes" << endl;
timeIntegrators[i].integrate();

Info<< "max(p): " << max(fluids[i].p()).value()
    << ", min(p): " << min(fluids[i].p()).value() << endl;
Info<< "max(T): " << max(fluids[i].T()).value()
    << ", min(T): " << min(fluids[i].T()).value() << endl;

//- Only update if viscosity is present
// if
// (
//     fluidRegions[i].foundObject<thermophysicalTransportModel>
//     (
//         thermophysicalTransportModel::typeName
//     )
// )
// {
//     devRhoReffs[i] ==
//         fluidRegions[i].lookupObject<thermophysicalTransportModel>(
//             thermophysicalTransportModel::typeName
//         ).momentumTransport().devTau();
// }

timeIntegrators[i].clear();
