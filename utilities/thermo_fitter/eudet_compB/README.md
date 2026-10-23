# EuDetCode output for Comp-B (no Al)

Products EOS and programmed-burn speed of the planar detonation tests
(`tests/1D/planar_znd`, `planar_znd_phaseTransfer`) and the reference of
their verifiers. Generated with the EuDetCode build in `external/eudetcode`
(2026-10-14):

```
CJDetonation "COMP-B" 1717 BKW -o eudet_compB
```

(1717 kg/m^3: just below the theoretical maximum density of the EuDetCode
COMP-B formulation, 1717.57 kg/m^3; 1720 is reset to the TMD.)
CJ state: D = 7934.30 m/s, p_CJ = 24.66 GPa, rho_CJ = 2224.5 kg/m^3,
T_CJ = 3891.4 K. The case entries follow from

```
python3 utilities/thermo_fitter/ejwlc_cj_entries.py utilities/thermo_fitter/eudet_compB
```

which sets the eJWLC E0 so that the products isentrope passes through the CJ
point for the case's reactant energy (Cv*300 K) and checks that the CJ
point of the resulting products EOS is at D (within 0.000%).
