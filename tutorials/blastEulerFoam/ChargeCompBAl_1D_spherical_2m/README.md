# 1D spherical ChargeCompBAl to 2 m (blastEulerFoam)

## Notes

The case of `tests/1D_spherical/ChargeCompBAl_1D_spherical_alCombustion`
(same `0`, `constant` and solver settings) on a 2 m domain, run to 2 ms, so
that the air blast is followed to Z ~ 2 m/kg^(1/3):

- Three phases: explosive `compBAl` (Murnaghan reactants, eJWLC products,
  EuDetCode products-only fit, D = 7154 m/s, p_CJ = 20.15 GPa), granular `Al`
  (30 um, 27% of the charge mass, burning with gas-phase oxidiser only) and
  the reacting `gas` (24 species, NASA-9 thermodynamics to 20 000 K,
  Abel-Nobel EOS, local equilibrium above 2500 K and a reversible GRI-3.0
  H2/O2/CO/C1 mechanism below).
- The products are transferred to the gas below 100 bar
  (`pressurePhaseChange`, tau = 20 us) and burn with the air.
- Charge radius 0.05 m (0.953 kg at 1820 kg/m3), initiated at the centre.
- Mesh: 200 uniform cells (0.5 mm, as in the test) to r = 0.1 m, then 599
  cells graded to 10 mm at r = 2 m. The outer boundary is not reached before
  2 ms.

Run with:
./Allrun

The run takes about 1 h 45 min on one core. `python3 summary.py` prints the
peak overpressure, arrival time and positive impulse at the pressure probes.

## Results (serial run, 2026-10-14)

| r [m] | peak overpressure [MPa] | arrival [us] |
|---|---|---|
| 0.08 | 42.5 | 11.1 |
| 0.10 | 32.7 | 14.3 |
| 0.15 | 20.5 | 23.8 |
| 0.30 | 7.93 | 65.1 |
| 0.50 | 3.39 | 148.5 |
| 1.00 | 0.759 | 511.5 |
| 1.50 | 0.266 | 1142.9 |
| 1.90 | 0.152 | 1812.9 |

On the graded mesh the shock is spread over several 5-10 mm cells beyond
1 m, so the peak lags the arrival (5% overpressure) by up to 90 us there and
is underestimated. Shocked air stays below 9300 K. The afterburning of the
products and the Al remains small in 1D: the products/air interface is
resolved but not mixed (multi-D mixing, Stage 5).
