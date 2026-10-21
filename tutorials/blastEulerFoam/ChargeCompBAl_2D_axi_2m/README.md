# 2D axisymmetric ChargeCompBAl to 2 m with AMR (blastEulerFoam)

## Notes

The full E2 setup of `ChargeCompBAl_1D_spherical_2m` (`setup_e2.py --variant
full`) on a 2D axisymmetric mesh, run to 2 ms (plan step G2):

- Comp-B/Al charge of radius 0.05 m (0.953 kg), initiated at the centre;
  explosive `compBAl` (eJWLC products), granular `Al` (30 um, eMelting, Loth
  drag, Kavanau heat transfer, Beckstead burning to Al2O3 with the gas
  oxidisers, ignition at 2050 K tracked by the transported ignited fraction)
  and the reacting `gas` (28 species, NASA-9, equilibrium above 2500 K,
  GRI-3.0 subset below). Products -> gas transfer below 100 bar (tau 20 us).
- Mesh: 5-degree wedge of the quarter plane x, y in [0, 2] m (axis y,
  symmetry plane y = 0, outer boundaries transmissive), base cells of 16 mm.
- Start from the 1D solution: the flow is spherically symmetric until the
  shock leaves the charge, and the detonation and breakout need 0.5 mm cells.
  `Allrun` runs the 1D case to 60 us in `oneD/` (skipped if `oneD/6e-05`
  exists), decomposes the 2D mesh and maps the 1D fields onto it with
  `rotateFields -refine -extend` in parallel (refining with
  `constant/dynamicMeshDict`). The mapping is done after `decomposePar`
  because `decomposePar` does not decompose the refinement history: cells
  refined before it can never be coarsened. `system/controlDict` sets
  `beginTime 0` so that the run counts as a restart (the products energy
  then does not get the detonation energy a second time).
- AMR (`hexRefiner`): refinement where neighbouring cells differ by more than
  20% in density or pressure (shock, products/air contact) and at the Al cloud
  edge (`scaledDelta`; `densityGradient` keeps the whole fireball refined).
  The finest level is a function of time and follows the cell size of the 1D
  graded mesh at the shock radius: 1 mm at the mapping, 2 mm to 350 us, 4 mm
  to 1.27 ms, then 8 mm (`maxLevel` table, floored; 0.5 mm before 40 us for
  a run from 0). Cells above the current level are coarsened.
- Probes on rays at 0, 45 and 90 degrees at the radii of the 1D case, and
  every 15 degrees at 0.3, 0.5, 1.0 and 1.5 m; sample lines at 0, 45 and 90
  degrees; peak pressure map (`radialPeakPressure`).

Allrun runs blastEulerFoam on 6 processes with load balancing (about 6 h;
18 000 - 25 000 cells, 6900 steps).

## Results (2026-10-21)

Against the 1D case at 0.3 - 1.9 m: arrival within 1%, peaks up to 9% and
positive impulse up to 27% (at 0.5 m) higher. The products mix with the air
through the fireball instabilities and afterburn more than in 1D (14% less CO
and 12% less H2 at 2 ms); the Al burns about as in 1D. Angular spread of the
peaks 1 - 5%, with the axis and the symmetry plane within the other
directions. The AMR uses 4% of the cell-steps of a uniform mesh at the same
finest levels.

Run with:
./Allrun
