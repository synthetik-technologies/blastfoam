fluidRegions[i].update();

fluids[i].encode();

Info<< "Calculating Fluxes" << endl;
timeIntegrators[i].integrate();

Info<< "max(p): " << max(fluids[i].p()).value()
    << ", min(p): " << min(fluids[i].p()).value() << endl;
Info<< "max(T): " << max(fluids[i].T()).value()
    << ", min(T): " << min(fluids[i].T()).value() << endl;

//- Only update if viscosity is present
if
(
    solidRegions[i].foundObject<blast::turbulenceModel>
    (
        blast::turbulenceModel::propertiesName
    )
)
{

    devRhoReffs[i] =
        solidRegions[i].lookupObject<blast::turbulenceModel>
        (
            blast::turbulenceModel::propertiesName
        ).devRhoReff();
}
