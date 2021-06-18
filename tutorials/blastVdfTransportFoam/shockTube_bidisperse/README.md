# Bidisperse shock tube using velocity distribution transport tutorial

## Notes

The OpenQBMM library is required to run this case.

This case uses two particle phases within a continuous air phase. This is a dilute simulation where packing of particles does not occur. The continuous gas phase is represented by air, and both particle phases have a density of 1300 kg/m^3 with a coefficient of restitution of 0.9. The diameters are 80e-6 m and 11e-6 m. The original case uses a small particle diamter of 10e-6, but this does not give correct results while 11e-6 m does, therefore this is assumed to be a typo in the original paper. Rather than the typical representation of particles using the two-fluid method, the transport of the size-velocity moments is used. We assume that each particle size has is own unique velocity and any velocity variance for a given size is negligible (mono-kinetic).

This method is not ideal for the given case, but it is meant to illustrate how to use the blastVdfTransportFoam solver.

Original plots were digitized using WebPlotDigitizer: https://github.com/ankitrohatgi/WebPlotDigitizer/releases The quality of the plots in the original publication (below) was not ideal for reproduction, but we have made an effort to reproduce them as faithfully as possible, and the points are found in the "validation/validationData" directory.

The calculation took approx. 20 s to run on a single core laptop at the reference resolution (e.g. 400 cells). Reference solution and plots from this run are in the "validation/referencePlots" directory. Plots will be automatically created using the postProcess and createGraphs utilities for subsequent runs.

## Reference

```
Lai, S., Houim, R.W., Oran, E.S., 2018. Effects of particle size and density on dust dispersion behind a moving shock. Phys. Rev. Fluids 3, 064306. https://doi.org/10.1103/PhysRevFluids.3.064306

```
