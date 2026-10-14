# 1D Spherical Stationary Charge Interface (blastEulerFoam)

## Overview
Reproducer for the unreacted charge/air interface of `ChargeCompBAl_1D_spherical`. The case is identical except that the charge is never initiated (`delay 1` on the detonation point) and it runs for 5 µs. The exact solution is the initial state: unreacted Comp-B/Al ($\alpha_{\text{compBAl}} = 0.818$, $\alpha_{\text{Al}} = 0.182$, $\rho = 1624.45$ / $2700\,\text{kg/m}^3$) at rest next to air at 1 atm.

It tests that the phase decode, the volume-fraction transport and the interphase coupling keep a stiff solid/gas material interface at rest without modifying the conserved variables (see `report/dilute_phase_conservation_plan.md`). Before the fix of the granular buoyancy term (mixture pressure instead of the individual fluid-phase pressures), this interface broke up within 1 µs with a conservative decode, sending trace explosive at GPa pressures and about 2000 m/s into the air.

## Execution
```bash
./Allrun
python3 ../verify_spherical_interface.py
```
