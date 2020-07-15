# Rate stick (JWL++) blastFoam Validation

## Notes

The presented set of cases are compared to Guilkey et al. (2007) and Souers et al. (2000) are used to validate the JWL++ model (unreacted: Murnaghan, reacted: JWL) by examining detonation velocity vs. inverse radius. Each case returns a calculated detonation velocity for a given rate stick diameter. The diameters of the rate sticks include 6 mm, 8 mm, 12.5 mm, 20 mm, and a one dimensional case (infinite radius). The rate sticks were activated using an initial velocity of 90 m/s at the bottom of the rate stick which created a high enough pressure to begin the activation process.

Original plots were digitized using WebPlotDigitizer: https://github.com/ankitrohatgi/WebPlotDigitizer/releases The quality of the plots in the original publication (below) was not ideal for reproduction, but we have made an effort to reproduce them as faithfully as possible, and the points are found in the "validation/validationData" directory.

All cases used a cell size of 0.25 mm. The run times on a four core desktop and number of cells for the cases were:

1D (400 cells): 17 s
6 mm (1920 cells): 90 s
8 mm (32000 cells): 3 min
12.5 mm (40000 cells): 4 min
20 mm (40000 cells): 4 min

Reference solution and plots from this run are in the "validation/referencePlots" directory. Plots will be automatically created using the postProcess and createGraphs utilities for subsequent runs. Results for the standard JWL and JWL-C model are provided in the reference plots, but the cases only use the JWL-C form.

## Reference

```
Guilkey, J.E., T.B. Harman, and B. Banerjee. “An Eulerian–Lagrangian Approach for Simulating Explosions of Energetic Devices.” Computers & Structures 85, no. 11–14 (June 2007): 660–74. https://doi.org/10.1016/j.compstruc.2007.01.031.
Souers, P.C., Anderson, S., Mercer, J., McGuire, E., Vitello, P., 2000. JWL++: A Simple Reactive Flow Code Package for Detonation 5.

```

