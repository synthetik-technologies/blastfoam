# 1D planar stationary charge (blastEulerFoam)

## Notes

`ChargeCompBAl_1D` on a coarser mesh (200 cells) with the detonation delayed
to t = 1 s, so the charge is never initiated within the 20 us run. The exact
solution is the initial state: unreacted Comp-B/Al (alpha.compBAl = 0.818,
alpha.Al = 0.182) at rest next to air at 1 atm. Any motion or pressure
change at the charge/air interfaces is a numerical error (see also
`tests/1D_spherical/ChargeCompBAl_1D_spherical_interface`).

Thermodynamics and phase coupling are those of `ChargeCompBAl_1D`.

Run with:
./Allrun
