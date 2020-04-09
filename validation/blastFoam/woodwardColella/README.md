# Woodward and Colella 1984 blastFoam Validation

## Notes

This problem was originally used by Woodward and Colella (1984) to compare the performance of several different hydrodynamical methods on problems involving strong shocks and narrow features. It has no analytical solution (except at very early times), but since it is one-dimensional, one may produce a converged solution by running the code with a large number of cells, permitting an estimate of the self-convergence rate.

Original plots were digitized using WebPlotDigitizer: https://github.com/ankitrohatgi/WebPlotDigitizer/releases The quality of the plots in the original publication (below) was not ideal for reproduction, but we have made an effort to reproduce them as faithfully as possible, and the points are found in the "validation/validationData" directory. The problem setup and grid size were as presented in the reference, and the blastFoam results match the published reference well.

The calculation took approx. 5 min to run on a single core laptop at the reference resolution (e.g. 9600 cells, 1m domain). Reference solution and plots from this run are in the "validation/referencePlots" directory. Plots will be automatically created using the postProcess and createGraphs utilities for subsequent runs.


## Reference

```
Woodward, Paul, and Phillip Colella. “The Numerical Simulation of Two-Dimensional Fluid Flow with Strong Shocks.” Journal of Computational Physics 54, no. 1 (1984): 115–173.
```

