# Under water detonation blastFoam validation

## Notes

This problem is meant to show the ability of the adaptive mesh refinement to correctly capture sharp density gradients using the error estimator presented in Zheng (2008). The case is initialized using a high pressure air bubble surrounded by water with a free surface.


The calculation took approx. 3 hours to run on a four core desktop at the reference resolution (e.g. 40x40 cells with a maximum refinement level of 5). Reference solution and plots from this run are in the "referencePlots" directory. The user run results must be viewed using a graphical viewer.


## Reference

```
Zheng, H.W., C. Shu, and Y.T. Chew. “An Object-Oriented and Quadrilateral-Mesh Based Solution Adaptive Algorithm for Compressible Multi-Fluid Flows.” Journal of Computational Physics 227, no. 14 (July 2008): 6895–6921. https://doi.org/10.1016/j.jcp.2008.03.037.

Shyue, K.-M., 1999. A Fluid-Mixture Type Algorithm for Compressible Multicomponent Flow with van der Waals Equation of State. Journal of Computational Physics 156, 43–88. https://doi.org/10.1006/jcph.1999.6349

```

