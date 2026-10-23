# 2D axisymmetric cylindrical ChargeCompBAl with AMR (blastEulerFoam)

## Notes

Plan step G3a (WP 5.2, `report/stage5_g3_plan.md`): a cylindrical Comp-B/Al
charge with the physics, mesh and refinement of `ChargeCompBAl_2D_axi_2m`
(the full E2 setup, 2 m wedge, AMR with the 2 mm mixing-layer criterion), to
compare the directional blast with the sphere of the same mass.

- Charge: cylinder about the axis (y) with the volume of the 0.05 m sphere of
  `ChargeCompBAl_2D_axi_2m` (0.953 kg), centred at the origin and initiated
  at its centre (programmed burn), so that the mid-plane y = 0 is a
  symmetry plane and the quarter plane x, y in [0, 2] m is kept. The 0
  degree ray is side-on, the 90 degree ray (the axis) end-on.
- `./Allrun <L/D>` (default 1) computes the radius and half-length, writes
  them into `system/setFieldsDict`, decomposes the mesh and sets and refines
  the charge to 0.5 mm with `setRefinedFields` in parallel (after
  `decomposePar`, so that the charge cells can be coarsened later; the
  uniform unburnt charge is coarsened in the first steps and refined again
  ahead of the detonation front). The charge volume on the mesh is 1.3%
  below the cylinder (staircase surface at L/D = 1).
- The programmed burn (`nodeProgrammed`) recomputes the cell arrival times
  after each topology change; arrival times mapped between refinement levels
  distorted the detonation front (35% ahead at 45 degrees at 2.5 us).
- The products/gas heat transfer is scaled by the reacted fraction
  (`field lambda.compBAl` in `constant/phaseProperties`): along the grazing
  detonation, hot shocked air at the charge surface otherwise heats traces of
  unreacted explosive to the temperature cap and the time step collapses.
- Explosive expanded below half its initial density is set fully reacted
  (`expansionBurnRatio 0.5` in `constant/physicalProperties.compBAl`):
  partly reacted explosive carried off the side by the grazing detonation
  otherwise expands to a few kg/m3, where the Murnaghan cold tension drives
  it to the temperature cap.
- Probes and sample lines as in `ChargeCompBAl_2D_axi_2m` (rays every 15
  degrees at 0.3, 0.5, 1.0 and 1.5 m).
- Mesh-dependence near the axis: the fireball within about 25 degrees of
  the axis is mesh-dependent (axis study in the report); the air shock on
  the axis agreed with 1D to a few per cent in the sphere run.

Runs on 6 processes with load balancing.

## Results

Run with:
./Allrun 1
