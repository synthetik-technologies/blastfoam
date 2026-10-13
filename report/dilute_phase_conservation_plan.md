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

- **Pressure and temperature relaxation in mixed cells (2026-10-13):** the inert and reacting spherical cases, the interface test and the 1D tutorials set `alphaDilute 1` (volume equilibration of every non-carrier fluid phase) and a `(compBAl in gas)` heat transfer of type `constantCoefficient` with `Kh 1e15`. `phaseSystem::relaxTemperature()` applies it as an exact exponential relaxation that conserves energy, so with this coefficient the temperatures equalise in one step. Result in the spherical reacting case: mixed cells within 1–15% in pressure and 10 K in temperature, and at 40 µs exactly in equilibrium (9.7 MPa, 3435 K). Planar tutorial at 20 µs: pressures within 1% and one common temperature, instead of products at 93 741 K beside gas at 124 K; the common temperature (about 9000 K) is still high. The interface test now measures the air criterion on cells without explosive or Al; the mixed interface cell (1.4e-6 explosive, warmed by 2 K, +1.4% pressure) is covered by the 10 bar bound.
  - `spherical_alCombustion` keeps `alphaDilute 0.01` and no products/gas temperature relaxation. Its Al burn releases mass into the gas phase inside the detonating charge, where the trace gas ($\alpha \sim 10^{-4}$) sits near its covolume limit (7000–8000 kg/m³, 1e10 Pa). Relaxing it to the products temperature crashes the gas chemistry at breakout (RKCK45 floating-point exception; with `seulex` the gas temperature inversion fails). `alphaDilute 1` alone also crashes, at 7.3 µs. Two further trials, both reverted: (i) no Al burning where the gas $\alpha < 10^{-3}$ (an oxidiser-phase gate in `reactingDiameterModel`) removes the in-charge burning that the verifier expects and still crashes at 2.5 µs, now from near-residual gas relaxed to the products temperature; (ii) skipping the gas chemistry where the gas $\alpha < 10^{-3}$ moves the crash to 7.1 µs. By then the gas inside the charge is pure Al vapour at $\alpha$ = 3–5e-3, 4500–5500 kg/m³ and 5 GPa, a state the Abel–Nobel/JANAF gas cannot represent. The underlying problem is the model: Al burnt inside the detonation products is released as vapour into the air phase. Resolved (2026-10-13) by burning only with gas-phase oxidiser: the `reactingDiameterModel` rate is scaled by the oxidiser-phase fraction of the surrounding fluid and its oxidiser mass fraction (`oxidiserPhase`, `oxidiserSpecies`, `YoxidiserRef`), so no Al burns where only products surround it. The case now runs to 40 µs with and without the p–T relaxation; it keeps the dilute-only setting because with the relaxation the unreacted charge/air interface becomes unstable at 2.9 µs (Section 6).

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
| `alphaDilute 1` (volume equilibration of every non-carrier fluid phase, 2026-10-13) | runs (spherical, interface, `ChargeCompBAl_1D`), pressures equal within 0.8–1.5 after 10 µs, probes almost unchanged; but at the products/air contact the carrier (largest $\alpha$) is the products phase, which is a trace by mass (0.24 vs 18 kg/m³ of air): it expands to 0.3 kg/m³ at 1e5 K while the air is compressed to 90 kg/m³ at 250 K |

## 5. Effect on the spherical charge results

`ChargeCompBAl_1D_spherical`, pressure probes (programmed burn, detonation velocity 7612 m/s):

| Probe | Original decode (`HEAD` before this work) | Conservative decode, without dilute equilibration | With dilute equilibration | Final (2026-10-13: no gas density floor, `alphaDeltaMin 1e-3`) | Products-only fit, α_Al = 0.182 (2026-10-13) |
|---|---|---|---|---|---|
| 1–4 cm (in the charge), peak | 11.1–13.5 GPa | 10.3–12.5 GPa | 10.3–12.5 GPa | 10.3–12.5 GPa | 13.7–17.0 GPa |
| 5 cm (charge edge), peak | 2.18 GPa | 5.54 GPa | 2.21 GPa | 1.45 GPa (1.40 GPa with p–T relaxation) | 2.85 GPa |
| 8 cm, arrival (p > 100 bar) | 16.2 µs | 10.5 µs | 10.7 µs | 12.5 µs | 11.3 µs |
| 8 cm, peak | 0.15 GPa (16.7 µs) | 0.49 GPa (17.2 µs) | 0.32 GPa (26.6 µs) | 33 MPa (30.6 µs, in the products) | 38 MPa (28.1 µs, in the products) |
| 10 cm, arrival | 25.4 µs | 15.0 µs | 19.4 µs | 17.9 µs (17.5 µs with p–T relaxation) | 15.0 µs |
| 10 cm, peak | 0.09 GPa (26.1 µs) | 0.27 GPa (31.7 µs) | 0.16 GPa (40 µs, end of run) | 14 MPa (18.4 µs, air shock) | 27 MPa (15.3 µs, air shock) |
| 15 cm, arrival / peak | – | – | not reached | 32.3 µs / 13 MPa | 25.7 µs / 19 MPa |
| Measured detonation velocity | 7650 m/s | 7640 m/s | 7640 m/s | 7661 m/s | 7664 m/s |

The air shock arrives much earlier without the resets: the original decode deleted the energy of the products expanding into the air. The averaged air-shock velocity over the first 3 cm is about 3 km/s with the resets and about 7 km/s without. With the gas density floor the 8 and 10 cm probes peaked well after the first arrival at 0.16–0.32 GPa: gas held at 1.2 kg/m³ while its conserved density was far lower. Without the floor the air shock is self-consistent: it runs from 10 to 15 cm at about 3.5 km/s, matching its 140 bar peak (M ≈ 11). At 8 cm the later maximum (33 MPa) occurs inside the expanded products. These near-field values have not been compared with experimental or reference data. The last column uses the products-only fit with the corrected Al fraction (before, the full fit with ρ0 = 1820 already contained the inert Al and the Al phase added α = 0.3 on top), with p–T relaxation.

## 6. Open

- **Non-dilute mixed cells:** fluid phases with $\alpha >$ `alphaDilute` still undergo the same volumetric strain (no pressure relaxation). In the spherical case (2026-10-13 state) about 50 cells at the products/air contact are mixed; there the products are at 2–5000 times the gas pressure, and the entrained air is rarefied (down to 4e-4 kg/m³) at its 5000 K ceiling. Equilibrating all phases (`alphaDilute 1`, Section 4) equalises the pressures but leaves the energy split between the phases wrong (products that are a trace by mass reach 1e5 K). Done for the inert and reacting cases with products/gas temperature relaxation (Section 3); not used in `spherical_alCombustion` (its interface was unstable with it; not retried since the 2026-10-14 fix), and the common temperature in the planar tutorial's mixed cells (about 9000 K) suggests numerical heating at the contact.
- **Interface start-up and instability (fixed 2026-10-14):** a phase that is a trace in a cell receives face pressure fluxes weighted by the face volume fraction of its dense neighbour, balanced by $P_I\nabla\alpha$ only if its face pressure equals the cell's interfacial pressure. In the HLLC flux the trace side had the dense phase's own impedance (a trace of explosive in the first air cell: 4.5e6 Pa·s/m), or vacuum where $\alpha \le 10^{-6}$ (the few-bar start-up wave). With `alphaTrace` (UPDATES.md Section 4.6) the trace side takes the mixture impedance: `_interface` stays below 1.4 mm/s and at 1.01 bar for 40 µs, and the detonating cases are unchanged except `_alCombustion`, whose interface wave is gone. Earlier fixes on the way: the instant velocity relaxation heated trace phases by their own velocity defect, and the granular Al energy lacked the compression work. Rejected: a ghost state (absent side = present side), which fails at the detonating contact, and a bracketed search for trace explosive in tension in `equilibrateDiluteVolumes()`.
