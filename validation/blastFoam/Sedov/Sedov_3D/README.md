# Sedov (3D) blastFoam validation

## Notes

This problem was introduced by Sedov (1959) to examine a point detonation. The case setup was taken from (Raga et al. 2012) and corresponds to a spherical point charge. A single cell is initialized with a pressure corresponding to a specified internal energy. This case is used to show the ability of the solver to handle strong shocks. The initial conditions are dependent on the mesh, and the amount of energy supplied to the system must be corrected if a different grid spacing is used.

Original plots were digitized using WebPlotDigitizer: https://github.com/ankitrohatgi/WebPlotDigitizer/releases The quality of the plots in the original publication (below) was not ideal for reproduction, but we have made an effort to reproduce them as faithfully as possible, and the points are found in the "validation/validationData" directory.. The problem setup and grid size were as presented in the reference, and the blastFoam results match the published reference well.

The calculation took approx. 10 min to run on a single core desktop at the reference resolution (e.g. 500 cells, 0.5m domain). Reference solution and plots from this run are in the "validation/referencePlots" directory. Plots will be automatically created using the postProcess and createGraphs utilities for subsequent runs.


## Reference

```
Raga, A. C., J. Cantó, L. F. Rodríguez, and P. F. Velázquez. “An Analytic Model for the Strong-/Weak-Shock Transition in a Spherical Blast Wave: Spherical Blast Wave.” Monthly Notices of the Royal Astronomical Society 424, no. 4 (August 21, 2012): 2522–27. https://doi.org/10.1111/j.1365-2966.2012.21208.x.


Sedov, L. I. “Similarity and Dimensional Methods in Mechanics.” Similarity and Dimensional Methods in Mechanics, New York: Academic Press, 1959, 1959. http://adsabs.harvard.edu/abs/1959sdmm.book.....S.

```

