# Tutorial for an FSI simulation of the effects of a detonation on a perpendicular flap in 2D

This tutorial is meant to show the use of the movingAdaptiveFvMesh class in combination with the preCICE interface and the OpenFOAM preCICE adapter. The solid mechanics solver CalculiX is used to solve the motion of the flap.

The case is set up on a 6 x 4 m domain with a 0.1 m flap at the center. A C4 charge with a radius of 0.05 is place 0.95 m from the flap at a height of 0.5 m from the ground.

case described in the [preCICE wiki](https://github.com/precice/precice/wiki/Tutorial-for-FSI-with-OpenFOAM-and-CalculiX).

It is known to work for CalculiX 2.13 with CGX, but it should also work with newer versions.

You may run the coupled simulation in serial using the script `Allrun` or in parallel with `Allrun -parallel`. The output of each step will be redirected to log files. You can cleanup the simulation using `Allclean`.

If you prefer to run the two simulations in two different terminals and watch their output on the screen, use the (simpler) scripts `runFluid` (or `runFluid -parallel`) and `runSolid`. Please always run the script `runFluid` first.

CGX is used to prepare the Solid participant and it should be executable by running `cgx`. If it has a different name (e.g. `cgx_2.13`), adapt the respective run script, or set an alias/link from `cgx` to `cgx_2.13`.

You may adjust the end time in the precice-config_*.xml, or interupt the execution earlier if you want.

The instructions to build the required preCICE libraries and adapters can be found below

[Get preCICE](https://github.com/precice/precice/wiki/Get-preCICE)

[Build the OpenFOAM adapter](https://github.com/precice/openfoam-adapter/wiki/Building)

[Build CalculiX](https://github.com/precice/calculix-adapter/wiki/Installation-instructions-for-CalculiX)

[Build Calculix adapter](https://github.com/precice/calculix-adapter/wiki/Building-the-Adapter)

## Disclaimer

This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM® and OpenCFD® trade marks.
