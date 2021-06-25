# Simple FSI flaps simulation

## Notes

This case is used to demonstrate the capability of blastFoam to handle fluid-structure interaction.

The blastFSIFoam solver is still under development and changing the case setup may result in crashes.

The case consists of a small c4 charge being detonated near a wall in 2 dimensions. The flap will deform when the blast wave reached it, and will effect how the wave moves past the flap. This is a modified version of the preCICE tutorial presented.

The calculation took approx. 2 hrs to run on a single core laptop at the reference resolution (e.g. 240x160 cells) and a maximum refinement level of 1.

