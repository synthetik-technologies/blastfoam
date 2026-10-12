# Conservative Treatment of Dilute Phases and Material Interfaces in blastEulerFoam

Status: **largely done** (2026-10-12; density floor and packing limit 2026-10-13). The decode resets are removed, and the conserved variables change only where a phase is residual. The 0D tests conserve energy to about 5e-9. The spherical verifiers, the new stationary-interface test and the blastEulerFoam tutorials pass. The open points are in Section 6.

## 1. Problem

`fluidPhaseModel::decode()` (`4cb49946`) and `granularPhaseModel::decode()` modified the conserved variables of dilute phases:

| Phase model | Action | Threshold |
|---|---|---|
| `fluidPhaseModel` | density blended towards `residualRho`, `alphaRho` rewritten | $\alpha < 10^{-3}$ |
| `fluidPhaseModel` | velocity blended to zero, `alphaRhoU` rewritten | $\alpha < 10^{-2}$ |
| `fluidPhaseModel` | energy reset to $h_e(288.15\,\text{K})$, `alphaRhoE` rewritten | $\alpha < 10^{-2}$ |
| `granularPhaseModel` | blended velocity/energy/$\Theta$; mass removed where $\alpha > \alpha_\text{max}$ | `residualAlpha`, `alphaMax` |

Mass, momentum and energy were not conserved. In `tests/0D/phaseTransfer_0D` and `_3phase`, the energy reset destroyed 1.8% and 1.1% of the total energy at 91 µs, when the products fell below 1% volume fraction. The same resets acted at every material interface and in every dilute region. In the spherical charge cases they deleted the energy of the products expanding into the air (Section 5).

## 2. Root causes found

With a fully conservative decode, the spherical Comp-B/Al cases produced spurious GPa pressure precursors ahead of the detonation and some runs crashed. The planar tutorial `ChargeCompBAl_1D` crashed at 18 µs. A stationary reproducer (Step 1) separated two causes.

1. **Granular buoyancy term.** The granular momentum equation used $-\alpha_s \sum_f \alpha_f \nabla p_f / (1-\alpha_s)$, built from the individual fluid-phase pressures. Upstream used $-\alpha_s\sum_f \nabla p_f$. A fluid phase's pressure is not defined where that phase is absent: the decode sets it to `small`. So at the unreacted charge/air interface, $\nabla p_\text{compBAl}$ and $\nabla p_\text{gas}$ were 1 atm jumps over one cell. The Al compacted and moved, and through $\sum\alpha = 1$ the fluid volume fractions were rescaled at constant $\alpha\rho$. That changed the decoded explosive density by percent levels, which the stiff Murnaghan EOS turns into GPa. With the charge made of compBAl only (no Al), the interface stayed at rest.
   **Fix:** $-\alpha_s \nabla p$ with the mixture fluid pressure $p$. Together with the fluid $p_I\nabla\alpha_k$ terms, this gives the mixture momentum balance $\partial_t(\rho\mathbf{u}) + \nabla p = \dots$ and preserves a uniform-pressure state.
2. **Trace fluid phases at fluid/fluid contacts.** The fluid $\alpha_k$ is advected with the mixture volumetric flux ($D\alpha_k/Dt = 0$) and there is no pressure relaxation, so all fluid phases in a cell undergo the same volumetric strain. In trace cells ($\alpha_k \sim 10^{-5}$–$10^{-2}$), $\alpha_k$ and $\alpha_k\rho_k$ also come from different discretisations. The decoded density $\alpha_k\rho_k/\alpha_k$ of a stiff phase then drifts: trace products ahead of the products/air contact reached 4000–10000 kg/m³ and 0.1–0.7 TPa.
   **Fix:** dilute fluid phases ($\alpha_k <$ `alphaDilute`, default 0.01) are brought towards pressure equilibrium with the carrier fluid phase by exchanging volume (Section 3).

## 3. Implemented (2026-10-12)

- **`fluidPhaseModel::decode()`:** conserved fields are modified only where $\alpha \le$ `residualAlpha` (reference state 288.15 K, as upstream). The velocity cap (20000 m/s) and the energy clamp to the thermo temperature range act on the primitives only.
- **`granularPhaseModel::decode()`:** the same treatment. Where $\alpha$ reaches `alphaMax`, the grain density becomes `alphaRho/alphaMax` instead of mass being removed. This is a last-resort bound: the packing limit must come from the frictional stress, which needs a small `alphaDeltaMin` (2026-10-13; see below).
- **Fluid density floor removed (2026-10-13):** `fluidPhaseModel::decode()` no longer applies `rho.max(residualRho)` outside residual cells. `residualRho` is the reference density of residual cells only. With the gas `residualRho 1.2`, the floor acted in 34–90 cells of the spherical runs (gas $\alpha$ up to 0.76, decoded density down to 0.06 kg/m³) and overstated their pressure up to twentyfold.
- **`alphaDeltaMin` (2026-10-13):** the JohnsonJackson frictional pressure is $F_r(\alpha-\alpha_\text{min})^\eta/\max(\alpha_\text{max}-\alpha, \alpha_{\Delta\min})^p$. With the upstream tutorial value `alphaDeltaMin 0.05` it is capped at about 10 kPa, so it cannot stop compaction and $\alpha$ runs into `alphaMax`. The original decode then deleted the excess mass (41% of the particles in `tests/1D/gunib1d_houim3`); the conservative decode stores it as overdense grains (Al up to 3620 kg/m³ at the centre of the spherical charge). `houim3` now uses $10^{-6}$ (compaction speed 183.8 m/s vs 185.25 m/s reference, mass conserved to 4e-12), and the Comp-B/Al tests and 1D tutorials use $10^{-3}$ (no cell reaches `alphaMax`; maximum Al density 3067 kg/m³).
- **`granularPhaseModel::solve()`:** buoyancy $-\alpha_s\nabla p$ from the mixture fluid pressure `fluid_.p()`.
- **`phaseSystem::equilibrateDiluteVolumes()`:** called in `decode()` after the fluid phases are decoded. For each cell and each dilute fluid phase $k$ (not the carrier $c$, the fluid phase with the largest $\alpha$):
  - A damped Newton iteration on the exchanged volume: $\Delta\alpha = (p_k - p_c)/(\rho_k c_k^2/\alpha_k + \rho_c c_c^2/\alpha_c)$, with backtracking while $|p_k - p_c|$ does not decrease.
  - Each increment moves the work $\bar p\,\Delta\alpha$ ($\bar p$ the mean of the two pressures) from one phase's energy to the other's.
  - Only $\alpha_k$, $\alpha_c$ and the two total energies change, so mass, momentum and total energy are conserved exactly. The phases are decoded again afterwards.
  - Bounds: the energies stay within the thermo temperature range. Per iteration the dilute phase may at most double its volume and the carrier may give at most 10% of its volume. Per call the dilute phase may at most halve its volume. Compressing trace gas inside a detonating explosive to tens of GPa takes it outside the Abel–Nobel/chemistry range: `spherical_alCombustion` crashed in the gas chemistry at breakout without this bound.
  - `alphaDilute` in `phaseProperties` (default 0.01; 0 disables).
- **New test `tests/1D_spherical/ChargeCompBAl_1D_spherical_interface`** (`verify_spherical_interface.py`, Tier 3): the spherical charge, never initiated, for 5 µs. Criteria: all phases at rest (< 5 m/s where $\alpha > 10^{-3}$), phase densities within 0.1%, $p < 10$ bar everywhere, air within 1% of $p_0$, mass and energy within 1e-5. Before fix 1 the interface broke up within 1 µs (velocities up to 2e4 m/s, compBAl density +590%, 2.8e8 Pa).

## 4. Approaches tried and rejected

| Approach | Result |
|---|---|
| Blending of primitives only (conserved fields untouched) | unstable within 2e-7 s: the fluxes are built from primitives inconsistent with the conserved state |
| Conservative pairwise relaxation of dilute-phase velocity and temperature | crashes; unexplained 2.5e-5 energy leak in `phaseTransfer_0D_3phase` |
| Dilute fluid phases take the carrier pressure (primitive only) | precursors persist, crashes: $p$ inconsistent with $\rho$ |
| Upstream `instant` pressure relaxation (`goodBroyden`) | $\Delta t$ collapses to 1e-15–1e-20 s (relaxed states with unphysical wave speeds), before and after fix 1 |
| Kapila compressibility term $\alpha_k(Z_W/Z_k - 1)\nabla\cdot\mathbf{u}$ in the fluid $\alpha$ equation | the planar tutorial crashed earlier (13 µs): the spurious $\nabla\cdot\mathbf{u}$ at the interface is absorbed by soft trace phases |
| Dilute-volume equilibration by bisection over a wide bracket | spurious roots (non-monotonic $p$ with the work term); energy grew 150-fold in the 0D test |
| Dilute-volume equilibration, one-sided (expansion only) | interface test failed: under-pressured trace gas in the charge reached the thermo clamps |

## 5. Effect on the spherical charge results

`ChargeCompBAl_1D_spherical`, pressure probes (programmed burn, detonation velocity 7612 m/s):

| Probe | Original decode (`HEAD` before this work) | Conservative decode, without dilute equilibration | With dilute equilibration | Final (2026-10-13: no gas density floor, `alphaDeltaMin 1e-3`) |
|---|---|---|---|---|
| 1–4 cm (in the charge), peak | 11.1–13.5 GPa | 10.3–12.5 GPa | 10.3–12.5 GPa | 10.3–12.5 GPa |
| 5 cm (charge edge), peak | 2.18 GPa | 5.54 GPa | 2.21 GPa | 1.45 GPa |
| 8 cm, arrival (p > 100 bar) | 16.2 µs | 10.5 µs | 10.7 µs | 12.5 µs |
| 8 cm, peak | 0.15 GPa (16.7 µs) | 0.49 GPa (17.2 µs) | 0.32 GPa (26.6 µs) | 33 MPa (30.6 µs, in the products) |
| 10 cm, arrival | 25.4 µs | 15.0 µs | 19.4 µs | 17.9 µs |
| 10 cm, peak | 0.09 GPa (26.1 µs) | 0.27 GPa (31.7 µs) | 0.16 GPa (40 µs, end of run) | 14 MPa (18.4 µs, air shock) |
| 15 cm, arrival / peak | – | – | not reached | 32.3 µs / 13 MPa |
| Measured detonation velocity | 7650 m/s | 7640 m/s | 7640 m/s | 7661 m/s |

The air shock arrives much earlier without the resets: the original decode deleted the energy of the products expanding into the air. The averaged air-shock velocity over the first 3 cm is about 3 km/s with the resets and about 7 km/s without. With the gas density floor the 8 and 10 cm probes peaked well after the first arrival at 0.16–0.32 GPa: gas held at 1.2 kg/m³ while its conserved density was far lower. Without the floor the air shock is self-consistent: it runs from 10 to 15 cm at about 3.5 km/s, matching its 140 bar peak (M ≈ 11). At 8 cm the later maximum (33 MPa) occurs inside the expanded products. These near-field values have not been compared with experimental or reference data.

## 6. Open

- **Non-dilute mixed cells:** fluid phases with $\alpha >$ `alphaDilute` still undergo the same volumetric strain (no pressure relaxation). At the products/air contact after breakout, products at $\alpha \sim 0.1$–0.4 reach the 100 K temperature floor at GPa pressures, and the gas reaches its `residualRho` floor and 5000 K ceiling. A robust full pressure relaxation remains to be done.
- **gunIB1D:** the local gunIB1D build stops with a floating-point exception on all three Houim cases, and `verify_gunib1d_benchmarks.py` points at `build/gunIB1D` and `Verification/`, which do not exist. The cross-solver comparison cannot be repeated.
- **Interface start-up:** the discrete interface sends a decaying acoustic wave of a few bar into the explosive (density change ~0.002%).
- `phaseTransfer_0D` WP 1.1 criterion $\Delta p/p < 0.01\%$ (currently 0.58–0.59%).
- Possible double counting of Al in the full Comp-B/Al cases (`UPDATES.md`, Section 7).
