# 2D Axisymmetric Spherical ChargeCompBAl (blastEulerFoam, AMR)

Plan step G1 (WP 5.1): the spherical charge of
`tests/1D_spherical/ChargeCompBAl_1D_spherical` on a 2D axisymmetric mesh,
to check that the multi-dimensional solver with adaptive mesh refinement
keeps the spherical symmetry and reproduces the 1D result.

## Setup
- Physics, initial state and end time (40 us) as in the 1D case: Comp-B/Al
  charge of radius 0.05 m (alpha.compBAl = 0.818, inert 30 um Al with
  alpha.Al = 0.182), origin-initiated programmed burn, JANAF / Abel-Nobel air.
- Mesh: 5-degree wedge about the y axis, quarter plane x, y in [0, 0.2] m
  with a symmetry plane at y = 0. Base cells of 2 mm, refined by `hexRefiner`
  by up to 2 levels to the 0.5 mm of the 1D case: the charge from the start
  (`setRefinedFields`), then the density gradient, pressure jumps and the Al
  cloud edge (`constant/dynamicMeshDict`).
- Probes on rays at 0, 30, 45, 60 and 90 degrees from the symmetry plane at
  the radii of the 1D probes; sample lines at 0, 45 and 90 degrees.
- Runs on 6 processes (`system/decomposeParDict`).

## Execution
```bash
./Allrun
python3 ../verify_2d_axi_spherical.py
```
The verifier needs the results of the 1D case (run its test first). It
checks the angular variation of the arrival time and peak pressure in the
air (30-60 degrees at 0.08 m, 10-80 degrees at 0.10 and 0.15 m: peak < 5%,
arrival < 2%) and compares the mean of the 30, 45 and 60 degree rays with
the 1D values (arrival < 1% in the air and 1.5% in the charge, peak < 3%).

The rays on the axis and the symmetry plane are reported, not checked: with
the HLLC flux the grid-aligned air shock has a carbuncle there (peaks up to
9% low on the axis, a trace of products running along both to the shock).
`hybridHLL yes` in `system/fvSchemes` (gas and compBAl) removes it but raises
the air-shock peaks by 2.5-3% at this resolution (UPDATES.md, Section 4.6).

Result (2026-10-20, 2.2 h on 6 processes): arrival -0.25 to -0.74% and peak
+0.9 / +2.6 / +2.7% against 1D at 0.08 / 0.10 / 0.15 m; angular variation
10-80 degrees: peak 0.3 / 3.0 / 2.5%, arrival 0.3 / 1.0 / 1.5%.
