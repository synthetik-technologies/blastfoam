# G3: Cylindrical charges and fireball mixing (Stage 5, WP 5.2–5.4)

Plan step G3 of `status_review_2026-10-14.md`: cylindrical charges (L/D) and
mixing (RMI/RTI, LES vs RANS, PaSR). Builds on G2/G2b
(`tutorials/blastEulerFoam/ChargeCompBAl_2D_axi_2m`: full E2 physics on a
2 m wedge with AMR and the 2 mm mixing-layer criterion) and on the axis study
of 2026-10-23 (the 0 and 90 degree rays are mesh-dependent in the fireball;
grid noise near the axis is a few per cent of the radial velocity).

## G3a: Cylindrical charges (WP 5.2)

- Case: `tutorials/blastEulerFoam/ChargeCompBAl_2D_axi_cylinder`, the G2
  tutorial with a cylindrical charge of the same volume (0.953 kg) about the
  axis, set and refined to 0.5 mm in parallel after `decomposePar`
  (`setRefinedFields`), initiated at the centre (the mid-plane is then a
  symmetry plane, so the quarter plane is kept). `./Allrun <L/D>`.
- Runs: L/D = 1 and 4 first (largest contrast), then 0.5 and 2; to 1.3 ms
  (the end of the G2b sphere run, which is the reference of the same mass and
  physics).
- Diagnostics: peak overpressure and positive impulse against the angle
  (every 15 degrees) at 0.3, 0.5, 1.0 and 1.5 m; side-on (0 degrees) over
  end-on (90 degrees) ratio; ratio to the sphere at the same radius; Al cloud
  and fireball shape.
- End initiation (half plane, twice the cells) after the centre-initiated
  set, if needed for the comparison with the literature.
- Criteria: the trends of the literature for centrally initiated cylinders
  (peak and impulse side-on rise with L/D and fall end-on; side-on peaks up
  to about three times the sphere's in the near field, approaching the
  sphere at large scaled distance) and, where data for comparable scaled
  distances exist, the side-on/end-on ratios. References: Knock & Davies,
  Shock Waves 23 (2013) (spherical equivalence of cylinders), Xiao, Andrae &
  Gebbeken (effect of charge shape and initiation, free air), Anastacio &
  Knock. Values to be taken from the papers before acceptance.
- End-on values are on the axis, where the fireball is mesh-dependent (axis
  study); the air shock there agreed with 1D to a few per cent in G2b, which
  is well below the expected side/end ratios.

## G3b: Seeded mixing-layer growth (WP 5.3)

- The G2b sphere with a seeded perturbation of the charge surface (single
  modes and multi-mode), amplitude well above the grid noise of the axis
  study, to measure the growth of the products/air mixing layer (RMI at
  breakout, RTI during the deceleration) and its effect on the CO and Al
  burnout; compare with mixing-layer growth reported for spherical and
  aluminised charges (e.g. Kuhl, Balakrishnan & Menon). Away from the axis
  (outside 25 degrees) only.

## G3c: Turbulence and combustion closures (WP 5.4)

- RANS (k-omega SST) and LES (WALE) on the G3b case, and PaSR for the gas
  chemistry; afterburning sensitivity to the closure.

## Status

- 2026-10-23: plan; G3a case set up.
