// Do not print out solver information
SolverPerformance<vector>::debug = 0;

solidModels[i].evolve();
solidModels[i].updateTotalFields();

// Turn solver information back on
SolverPerformance<vector>::debug = 1;

//- Clear global Patches since displacement may have changed
globalPolyBoundaryMesh::New(solidRegions[i]).movePoints();

