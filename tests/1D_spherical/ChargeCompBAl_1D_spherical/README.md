# 1D Spherical ChargeCompBAl (blastEulerFoam)

## Overview
This test case models a 1D spherical sector detonation and multiphase blast wave expansion initiated at the origin ($r = 0$). It is based on `tutorials/blastEulerFoam/ChargeCompBAl_1D`, formulated on an axisymmetric wedge mesh radiating from the origin, and calibrated to the `EuDetCode` ideal detonation state for aluminized Composition B (70% Comp-B / 30% Al).

## Key Parameters
- **Solver**: `blastEulerFoam` (three-phase Eulerian system: reacting explosive `compBAl`, granular particulate `Al`, and surrounding ambient `gas`).
- **Domain**: 1D spherical wedge radiating from $(0, 0, 0)$ along the $X$-axis with outer radius $R = 0.2\,\text{m}$.
- **Charge Geometry**: Spherical explosive charge of radius $R_{\text{charge}} = 0.05\,\text{m}$ located at the origin ($r \in [0, 0.05]\,\text{m}$), initialized with $\alpha_{\text{compBAl}} = 0.7$, $\alpha_{\text{Al}} = 0.3$, $\alpha_{\text{gas}} = 0$.
- **Initiation**: Origin-initiated programmed detonation with `points ((0 0 0))` at $t = 0$.
- **Thermodynamics (`EuDetCode` Calibration)**:
  - Initial density: $\rho_0 = 1820\,\text{kg/m}^3$
  - Chapman-Jouguet Detonation Velocity: $D_{\text{CJ}} = 7612.3\,\text{m/s}$
  - Chapman-Jouguet Detonation Pressure: $P_{\text{CJ}} = 2.153 \times 10^{10}\,\text{Pa}$ ($21.53\,\text{GPa}$)
  - JWL Parameters: $A = 1732.51\,\text{GPa}$, $B = 11.717\,\text{GPa}$, $C = 0.8317\,\text{GPa}$, $R_1 = 5.8162$, $R_2 = 1.51614$, $\omega = 0.24273$, $C_v = 1410.19\,\text{J/(kg K)}$, $E_0 = -3.651 \times 10^6\,\text{J/kg}$.

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
