# Axisymmetric ChargeCompBAl (blastEulerFoam)

## Notes

A cylindrical Comp-B/Al charge (radius 0.26 m, height 0.26 m) standing on
the ground, initiated at the centre of its base and surrounded by
atmospheric air. The mesh is a 5-degree wedge, 1 m x 1 m, with 2.5 mm cells.

Thermodynamics and phase coupling are those of
`tests/1D_spherical/ChargeCompBAl_1D_spherical`:

- Three phases: explosive `compBAl` (Murnaghan reactants, eJWLC products),
  inert granular `Al` (1 mm particles) and the surrounding `gas` (air,
  76% N2 / 24% O2 by mass, as a multi-species JANAF / Abel-Nobel gas).
- The explosive uses the EuDetCode products-only eJWLC fit
  (`utilities/thermo_fitter/products_only_fit.py`): rho0 = 1624.45 kg/m3,
  D = 7154.23 m/s, p_CJ = 20.15 GPa. The 27% (by mass) Al is the separate
  `Al` phase: alpha.compBAl = 0.818, alpha.Al = 0.182 (charge density
  1820 kg/m3).
- Mixed products/air cells are kept in pressure and temperature equilibrium
  (`alphaDilute 1` and a stiff `compBAl in gas` heat transfer); trace phases
  use the mixture impedance at faces (`alphaTrace 1e-3`).

Allrun runs blastEulerFoam in parallel on 4 processes
(`system/decomposeParDict`) to 200 us.

Run with:
./Allrun
