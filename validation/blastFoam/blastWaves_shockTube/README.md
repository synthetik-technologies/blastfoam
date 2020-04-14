# Woodward and Colella 1984 blastFoam Validation

## Notes

This simple shock tube problem was presented by Needham (2018) and compared to the results from the SHAMRC code ("SHAMRC Capabilities and Description").

Original plots were digitized using WebPlotDigitizer: https://github.com/ankitrohatgi/WebPlotDigitizer/releases The quality of the plots in the original publication (below) was not ideal for reproduction, but we have made an effort to reproduce them as faithfully as possible, and the points are found in the "validation/validationData" directory.

The calculation took approx. 5 sec to run on a single core laptop at the reference resolution (e.g. 2000 cells, 250m domain). Reference solution and plots from this run are in the "validation/referencePlots" directory. Plots will be automatically created using the postProcess and createGraphs utilities for subsequent runs.


## Reference

```
Needham, Charles E. Blast Waves. Shock Wave and High Pressure Phenomena. Cham: Springer International Publishing, 2018. https://doi.org/10.1007/978-3-319-65382-2.

Needham, Charles, Joseph Crepeau, and Hank Happ. “SHAMRC Capabilities and Description.” n.d.

```

