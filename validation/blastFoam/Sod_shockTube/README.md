# Sod shock tube blastFoam Validation

## Notes

This case is a the standard Sod shock tube initially presented by Sod (1978).

Original plots were digitized using WebPlotDigitizer: https://github.com/ankitrohatgi/WebPlotDigitizer/releases The quality of the plots in the original publication (below) was not ideal for reproduction, but we have made an effort to reproduce them as faithfully as possible, and the points are found in the "validation/validationData" directory.

The calculation took approx. 0.5 s to run on a single core laptop at the reference resolution (e.g. 500 cells). Reference solution and plots from this run are in the "validation/referencePlots" directory. Plots will be automatically created using the postProcess and createGraphs utilities for subsequent runs.


## Reference

```
Sod, Gary A. “A Survey of Several Finite Difference Methods for Systems of Nonlinear Hyperbolic Conservation Laws.” Journal of Computational Physics 27, no. 1 (April 1, 1978): 1–31. https://doi.org/10.1016/0021-9991(78)90023-2.

```

