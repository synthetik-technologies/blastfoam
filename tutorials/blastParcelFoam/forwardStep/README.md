# blastParcelFoam forwardStep tutorial

## Notes

This is a modified version of the standard forwardStep tutorial from OpenFOAM. The case has been modified to include inertial, heated, Lagrangian particles.

Two methods to view the lagrangian particles exist: The standard OpenFOAM paraview reader can be used (*.OpenFOAM extension), but requires the case to first be reconstructed using reconstructPar. The particles can the be viewed using the "Glyph" filter. The second method is to use the "convertLagrangianPositions" to modifiy the outputs of the solver so that the standard paraview OpenFOAM (*.foam extension) can read the lagrangian files without. If you would like to restart the simulation after using this method, you must use the "-revert" flag so that OpenFOAM can correclty read in the lagrangian fields. The second method is used for this case.

The case takes approximately 30 min to run on a four-core desktop.


