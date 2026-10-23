# 1D Spherical ChargeCompBAl (blastEulerFoam)

## Overview
This test case models a 1D spherical sector detonation and multiphase blast wave expansion initiated at the origin ($r = 0$). It is based on `tutorials/blastEulerFoam/ChargeCompBAl_1D`, formulated on an axisymmetric wedge mesh radiating from the origin, and calibrated to the `EuDetCode` ideal detonation state for aluminized Composition B (Comp-B with 27% Al by mass).

## Key Parameters
- **Solver**: `blastEulerFoam` (three-phase Eulerian system: reacting explosive `compBAl`, granular particulate `Al`, and surrounding ambient `gas`).
- **Domain**: 1D spherical wedge radiating from $(0, 0, 0)$ along the $X$-axis with outer radius $R = 0.2\,\text{m}$.
- **Charge Geometry**: Spherical explosive charge of radius $R_{\text{charge}} = 0.05\,\text{m}$ located at the origin ($r \in [0, 0.05]\,\text{m}$), initialized with $\alpha_{\text{compBAl}} = 0.818$, $\alpha_{\text{Al}} = 0.182$, $\alpha_{\text{gas}} = 0$.
- **Initiation**: Origin-initiated programmed detonation with `points ((0 0 0))` at $t = 0$.
- **Thermodynamics (`EuDetCode` Calibration)**:
  - Products-only eJWLC fit (EuDetCode `--jwl-products-only`, `utilities/thermo_fitter/products_only_fit.py`): the 27% (by mass) inert Al is carried by the `Al` phase, so the explosive phase holds the active products only.
  - Reference density of the active phase: $\rho_0 = 1624.45\,\text{kg/m}^3$ (charge density $1820\,\text{kg/m}^3$ with $\rho_{\text{Al}} = 2700\,\text{kg/m}^3$)
  - Chapman-Jouguet Detonation Velocity: $D_{\text{CJ}} = 7154.23\,\text{m/s}$
  - Chapman-Jouguet Detonation Pressure: $P_{\text{CJ}} = 2.015 \times 10^{10}\,\text{Pa}$ ($20.15\,\text{GPa}$)
  - eJWLC Parameters: $A = 224.943\,\text{GPa}$, $B = 2.77303\,\text{GPa}$, $C = 1.47024\,\text{GPa}$, $R_1 = 3.456$, $R_2 = 0.813403$, $\omega = 0.302717$, $C_v = 1450\,\text{J/(kg K)}$, $T_{\text{CJ}} = 3844.39\,\text{K}$, $\rho_{\text{CJ}} = 2150.55\,\text{kg/m}^3$, $\omega_T = 0.371065$, $E_0 = -5.0347 \times 10^5\,\text{J/kg}$ (on the JANAF sensible-energy datum of the gas phase); programmed-burn $e_0 = 5.0174 \times 10^6\,\text{J/kg}$.
  - Gas phase: air (76% N2, 24% O2 by mass) as a multi-species JANAF / Abel-Nobel gas carrying the detonation-product species.

## Execution
Run the case directly using:
```bash
./Allrun
```
Clean the case with:
```bash
./Allclean
```
Verify physical accuracy and CJ front propagation speed:
```bash
python3 ../verify_spherical_compBAl.py
```
