# movingCone tutorial

## Notes

This case is a modified version of the OpenFOAM movingCone tutorial. The main modification is the use of the movingAdaptiveFvMesh dynamic mesh class. This allows for automatic mesh refinement in addition to the mesh movement.

If you would like to run this case in parallel, the Synthetik branch of OpenFOAM must be used due to a bug that causes moving, adaptive meshes to fail in parallel.

The case takes approximately 1 hour to run on a single core desktop.
