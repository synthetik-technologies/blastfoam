# Validation of polydispserse granular

## Notes

This case uses six granular phases within a continuous air phase. This a dilute simulation where packing of particles does no occur. The continuous gas phase is represented by air. Five of the particles phases are used to represent a continuous size distribution with a density 2700 kg/m^3 and a coefficient of restitution of 0.9. The last granular phase is that of coal with a density of 1300 kg/m^3, coefficient of restitution of 0.9, and a diameter of 10e-6 m.

Original plots were digitized using WebPlotDigitizer: https://github.com/ankitrohatgi/WebPlotDigitizer/releases The quality of the plots in the original publication (below) was not ideal for reproduction, but we have made an effort to reproduce them as faithfully as possible, and the points are found in the "validation/validationData" directory.

The calculation took approx. 3 min to run on a single core laptop at the reference resolution (e.g. 500 cells). Reference solution and plots from this run are in the "validation/referencePlots" directory. Plots will be automatically created using the postProcess and createGraphs utilities for subsequent runs.

## Reference

```
Lai, S., Houim, R.W., Oran, E.S., 2018. Effects of particle size and density on dust dispersion behind a moving shock. Phys. Rev. Fluids 3, 064306. https://doi.org/10.1103/PhysRevFluids.3.064306


```
