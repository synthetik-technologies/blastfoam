solidModels[i].evolve();

//- Clear global Patches since displacement may have changed
globalPolyBoundaryMesh::New(solidRegions[i]).movePoints();

