# 2-D Riemann problem blastFoam validation

## Notes

This problem was simulated by Lax and Liu (1998) The interaction of multiple shocks occuring in two dimensions. The case consists of four quadrants initialize with varying densities, pressures, and velocities.

Original plots were digitized using WebPlotDigitizer: https://github.com/ankitrohatgi/WebPlotDigitizer/releases The quality of the plots in the original publication (below) was not ideal for reproduction, but we have made an effort to reproduce them as faithfully as possible, and the points are found in the "validation/validationData" directory.. The problem setup and grid size were as presented in the reference, and the blastFoam results match the published reference well.

The calculation took approx. 6 min to run on a four core desktop at the reference resolution (50 x 50 cells) on a 1 m  x 1 m domain with adaptive mesh refinement with a mamximum refinement level of 3. A reference solution from this run are in the "validation/referencePlots" directory. The user run results must be viewed using a graphical viewer.


## Reference

```
Lax, Peter D., and Xu-Dong Liu. “Solution of Two-Dimensional Riemann Problems of Gas Dynamics by Positive Schemes.” SIAM Journal on Scientific Computing 19, no. 2 (March 1998): 319–40. https://doi.org/10.1137/S1064827595291819.
```

