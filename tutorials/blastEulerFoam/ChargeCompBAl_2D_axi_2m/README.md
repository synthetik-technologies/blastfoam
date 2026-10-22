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
  The products/air mixing layer is held at 2 mm (`mixingLevel`) by the
  `linearCombination` criterion on the carbon/(carbon + nitrogen) element
  ratio of the gas (0.449 in the products, 0 in air; refined between 0.01 and
  0.44): once the layer is spread over a few cells, the jumps of density and
  pressure no longer mark it.
- Probes on rays at 0, 45 and 90 degrees at the radii of the 1D case, and
  every 15 degrees at 0.3, 0.5, 1.0 and 1.5 m; sample lines at 0, 45 and 90
  degrees; peak pressure map (`radialPeakPressure`).

Allrun runs blastEulerFoam on 6 processes with load balancing (about 6 h;
18 000 - 25 000 cells, 6900 steps).

## Results

2026-10-21, without the mixing-layer criterion (to 2 ms, 5.6 h): the mixing
layer diffused on the 16 mm base cells (3 - 4 times wider than in 1D), the
products afterburnt more than in 1D (14% less CO at 2 ms) and peaks and
positive impulse were up to 9% and 27% higher.

2026-10-22, with the mixing layer at 2 mm (to 1.3 ms, 18 h on 6 processes,
up to 71 000 cells): CO and H2 burnout within 1% of 1D, peaks within -3 /
+4% and positive impulse within -4 / +7% of 1D at 0.3 - 1.5 m, arrival
within 2.5%. From 0.5 ms the fireball surface forms Rayleigh-Taylor fingers
(about 50 - 60 mm apart), seeded by the mesh, so not mesh-converged. Along
the axis air runs into the fireball and along the symmetry plane a jet of
products runs out (inward secondary shock ahead along both; under
investigation). A run to 2 ms takes about 40 h.

Run with:
./Allrun
