const fvMesh& mesh = solidRegions[i];

solidThermoModel& thermo = solidThermos[i];

volVectorField& D = Ds[i];
pointVectorField& pointD = pointDs[i];
volVectorField& divSigmaExp = divSigmaExps[i];
volSymmTensorField& sigmaD = sigmaDs[i];

// Mechanical fields
const volScalarField& mu = mus[i];
const volScalarField& nu = nus[i];
const volScalarField& lambda = lambdas[i];
const volScalarField& E = Es[i];
Switch& compactNormalStress = compactNormalStresses[i];
const Switch& planeStress = planeStresses[i];


// Thermal fields
Switch& thermalStress = thermalStresses[i];


fv::options& fvOptions = solidSources[i];
