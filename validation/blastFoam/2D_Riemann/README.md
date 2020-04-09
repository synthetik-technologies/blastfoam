# 2-D Riemann problem blastFoam validation

## Notes

This problem was simulated by Lax and Liu (1998) The interaction of multiple shocks occurring in two dimensions. The case consists of four quadrants initialize with varying densities, pressures, and velocities.

The calculation took approx. 6 min to run on a four core desktop at the reference resolution (50 x 50 cells) on a 1 m  x 1 m domain with adaptive mesh refinement with a maximum refinement level of 3. A reference solution from this run are in the "validation/referencePlots" directory. The user run results must be viewed using a graphical viewer.


## Reference

```
Lax, Peter D., and Xu-Dong Liu. “Solution of Two-Dimensional Riemann Problems of Gas Dynamics by Positive Schemes.” SIAM Journal on Scientific Computing 19, no. 2 (March 1998): 319–40. https://doi.org/10.1137/S1064827595291819.
```

