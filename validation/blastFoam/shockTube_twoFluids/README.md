# Two fluid shock tube detonation blastFoam validation

## Notes

This simple shock tube problem is used to validation the two-fluid solver consisting of fluid phase phase described by a stiffened equation of state and a gas phase described by the Van der Walls equation of state. The case is compared to the results of Zheng (2011).

Original plots were digitized using WebPlotDigitizer: https://github.com/ankitrohatgi/WebPlotDigitizer/releases The quality of the plots in the original publication (below) was not ideal for reproduction, but we have made an effort to reproduce them as faithfully as possible, and the points are found in the "validation/validationData" directory.

The calculation took approx. 1 sec to run on a single core desktop at the reference resolution (e.g. 300 cells). Reference solution and plots from this run are in the "referencePlots" directory. Plots will be automatically created using the postProcess and createGraphs utilities for subsequent runs.


## Reference

```
Zheng, H. W., C. Shu, Y. T. Chew, and N. Qin. “A Solution Adaptive Simulation of Compressible Multi-Fluid Flows with General Equation of State.” International Journal for Numerical Methods in Fluids 67, no. 5 (2011): 616–37. https://doi.org/10.1002/fld.2380.

```

