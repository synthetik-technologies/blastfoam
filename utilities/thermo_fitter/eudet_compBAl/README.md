# EuDetCode output for 70% Comp-B / 30% Al

Inputs of the Comp-B/Al fits and of the report figures (`utilities/generate_report_figures.py`).
Generated with the EuDetCode build in `external/eudetcode` (2026-10-14):

```
CJDetonation "RDX:0.42 TNT:0.28 Al:0.3" 1820 BKW --inert "Al(s):0.9" --jwl-products-only -o products
CJDetonation "RDX:0.42 TNT:0.28 Al:0.3" 1820 BKW --inert "Al(s):0.9" -o full
```

`products/` is the fit used in the cases (D = 7154.23 m/s, p_CJ = 20.15 GPa);
`full/` is the full-mixture fit, for comparison only. The case entries follow from

```
python3 utilities/thermo_fitter/products_only_fit.py utilities/thermo_fitter/eudet_compBAl/products tests/0D/phaseTransfer_0D
```
