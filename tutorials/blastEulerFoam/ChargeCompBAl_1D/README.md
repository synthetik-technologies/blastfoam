# 1D planar ChargeCompBAl (blastEulerFoam)

## Notes

Planar 1D version of the Comp-B/Al charge, with the same thermodynamics and
phase coupling as `tests/1D_spherical/ChargeCompBAl_1D_spherical`:

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

The domain is 0.2 m long along the X-axis (1000 cells). The charge fills
X = 0.01-0.15 m and is initiated at X = 0.0505 m. The left boundary is
pressureWaveTransmissive, the right boundary zeroGradient. The run ends at
20 us.

Run with:
./Allrun
