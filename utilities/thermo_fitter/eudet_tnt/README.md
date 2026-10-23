# EuDetCode output for TNT

Products EOS and programmed-burn speed of the TNT spherical blast test
(`tests/1D_spherical/tnt_blast`, comparison with Kingery-Bulmash).
Generated with the EuDetCode build in `external/eudetcode` (2026-10-15,
run from `external/eudetcode`, which holds the `Data/` directory):

```
CJDetonation "TNT" 1630 BKW -o <this directory>
```

CJ state: D = 6696.23 m/s, p_CJ = 16.89 GPa, rho_CJ = 2119.74 kg/m^3,
T_CJ = 3636.8 K (BKW; measured for cast TNT at 1.63 g/cm^3 about 6.9 km/s
and 19-21 GPa). The case entries follow from

```
python3 utilities/thermo_fitter/ejwlc_cj_entries.py utilities/thermo_fitter/eudet_tnt --e-reactants 300000
```

with the reactant energy Cv*300 K for the EuDetCode TNT heat capacity
(227.1 J/(mol K) = 1000 J/(kg K)); E0 = -4262687.88 J/kg. The unreacted TNT
is Murnaghan with the EuDetCode ingredient data K0 = 8.93 GPa, n = 5.5.
