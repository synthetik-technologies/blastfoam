# Sedov (2D) blastFoam validation

## Notes

This problem was introduced by Sedov (1959) to examine a point detonation. The case setup was taken from ("34.1 Hydrodynamics Test Problems") and corresponds to a cylinder. A single cell is initialized with a pressure corresponding to a specified internal energy. This case is used to show the ability of the solver to handle strong shocks. The initial conditions are dependent on the mesh, and the amount of energy supplied to the system must be corrected if a different grid spacing is used.

Original plots were digitized using WebPlotDigitizer: https://github.com/ankitrohatgi/WebPlotDigitizer/releases The quality of the plots in the original publication (below) was not ideal for reproduction, but we have made an effort to reproduce them as faithfully as possible, and the points are found in the "validation/validationData" directory.. The problem setup and grid size were as presented in the reference, and the blastFoam results match the published reference well.

The calculation took approx. 40 s to run on a single core desktop at the reference resolution (e.g. 400 cells, 0.4m domain). Reference solution and plots from this run are in the "validation/referencePlots" directory. Plots will be automatically created using the postProcess and createGraphs utilities for subsequent runs.


## Reference

```
“34.1 Hydrodynamics Test Problems.” Accessed April 8, 2020. http://flash.uchicago.edu/site/flashcode/user_support/flash_ug_devel/node184.html#SECTION010114000000000000000.

Sedov, L. I. “Similarity and Dimensional Methods in Mechanics.” Similarity and Dimensional Methods in Mechanics, New York: Academic Press, 1959, 1959. http://adsabs.harvard.edu/abs/1959sdmm.book.....S.

```

