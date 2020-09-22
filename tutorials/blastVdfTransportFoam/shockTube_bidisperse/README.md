# Validation of bidisperse granular

## Notes

This case uses two granular phases within a continuous air phase. This a dilute simulation where packing of particles does no occur. The continuous gas phase is represented by air, and the both particle phases have a density of 1300 kg/m^3 with a coefficient of restitution of 0.9. The diameters are 80e-6 m and 11e-6 m. The original case uses a small particle diamter of 10e-6, but this does not give correct results while 11e-6 m does, therefore this is assumed to be a typo in the original paper.

In contrast to the original case, no velocity variance in considered for this case resulting in less diffusion of particles. Additionally, the velocity outside of the areas with particles are non-zero due to the way in which the mean velocities are computed.

Original plots were digitized using WebPlotDigitizer: https://github.com/ankitrohatgi/WebPlotDigitizer/releases The quality of the plots in the original publication (below) was not ideal for reproduction, but we have made an effort to reproduce them as faithfully as possible, and the points are found in the "validation/validationData" directory.

The calculation took approx. 20 s to run on a single core laptop at the reference resolution (e.g. 400 cells). Reference solution and plots from this run are in the "validation/referencePlots" directory. Plots will be automatically created using the postProcess and createGraphs utilities for subsequent runs.

## Reference

```
Lai, S., Houim, R.W., Oran, E.S., 2018. Effects of particle size and density on dust dispersion behind a moving shock. Phys. Rev. Fluids 3, 064306. https://doi.org/10.1103/PhysRevFluids.3.064306


```
