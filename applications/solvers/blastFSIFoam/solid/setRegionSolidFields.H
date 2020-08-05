const fvMesh& mesh = solidRegions[i];

solidThermoModel& thermo = solidThermos[i];

volVectorField& D = Ds[i];
volVectorField& divSigmaExp = divSigmaExps[i];
volSymmTensorField& sigmaD = sigmaDs[i];

// Mechanical fields
volScalarField& mu = mus[i];
volScalarField& lambda = lambdas[i];
Switch& compactNormalStress = compactNormalStresses[i];


// Thermal fields

volScalarField& threeKalpha = threeKalphas[i];
volScalarField& DT = DTs[i];
Switch& thermalStress = thermalStresses[i];


fv::options& fvOptions = solidSources[i];
