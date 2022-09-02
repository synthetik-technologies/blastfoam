# mapped building 3D

## Notes

These cases are to show how to use the rotateFields utility to run a refined 1D/2D case to initialize a refined solution and then rotate to a 3D case. The sector (1D, serial) initialization is used by default, but the wedge (2D, parallel) initializtion can be used by running './Allwmake -wedge'. Both methods of initialization will write to the symbolicly linked 'axisymmetricCharge' folder.

The sector set of cases takes approximately 2 minutes to run the wedge case take approximately 10 minutes, both on a four core desktop.
