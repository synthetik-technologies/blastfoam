# Status Review and Forward Plan (2026-10-14)

Review of the code, tests and documents at `ab275c9e` (branch `develop`) against the aims of `report/TBX-blastFoam.tex` and the staged plan in `report/development_plan.md`. Every statement below was checked against the code or the test cases, not only against the documents.

---

## 1. The aim

The report (Section 1) sets out to simulate an aluminised explosive (Comp-B/Al, 70/30) from initiation to late-time afterburning in air, covering four phases:

1. detonation and early product expansion with a high-pressure EOS;
2. mixing of the products with air, with finite-rate chemistry in the mixture;
3. transport, heating, ignition and combustion of the Al particles;
4. afterburning of Al in the mixture of products and air.

The end point (Stage 6, Gate 5) is a 3D `ChargeCompBAl` simulation whose peak overpressure, arrival time and positive impulse lie within 8% of experimental gauge data. A downselected dispersed-phase method and production guidelines go with it.

---

## 2. Overall status

Phase (1) works. Phases (2)–(4) are not yet represented in the solver runs: the code exists for parts of them, but no test case exercises them together.

| Report phase | What exists and is verified | What is missing |
|---|---|---|
| (1) Detonation and early expansion | Programmed burn with eJWLC products-only fit (EuDetCode, D = 7154 m/s); 1D spherical 3-phase charge in blastEulerFoam, D within 0.12%; conservative multi-fluid decode; stable charge/air interface (`alphaTrace`); constant-$C_p$ high-T air | Peak in charge 16.6 GPa vs $p_\text{CJ}$ 20.15 GPa (resolution, 400 cells); air above 5000 K is too hot (no dissociation) |
| (2) Products/air mixing and chemistry | `pressurePhaseChange` products → reactive gas, verified in 0D (energy 4e-9, $\Delta p/p$ 3e-7) and in planar ZND; gas-phase chemistry (CH4, CO, H2, Al reactions) | **Phase transfer is disabled in every Comp-B/Al spherical case and tutorial** (`constant/fvModels.disabled`). The products stay in the eJWLC phase, the gas contains only air (`CO.gas = CH4.gas = 0`), so no product afterburning can occur. |
| (3) Al heating, ignition, combustion | Granular Al phase with drag, heat transfer and compaction; `reactingParticle` mass transfer with the empirical `pressureBased` rate, restricted to gas-phase oxidiser | Melting, ignition threshold, Beckstead $d^n$ multi-oxidiser rate, Mach-dependent drag and alumina condensation exist **only as numpy models** (`tests/0D/verify_al_particle_kinetics.py`), not in the solver. Al reaches 578 K by 40 µs, so nothing ignites. |
| (4) Afterburning in air | Al + O2/CO2/H2O gas reactions in `_alCombustion` | 0.0044% of the Al burns by 40 µs; the domain is 0.2 m and 40 µs, while afterburning takes metres and milliseconds |

### 2.1 Stage-by-stage

| Stage | Plan status | Assessment |
|---|---|---|
| 0 Toolchain | Partial | Harness (18 tests, 4 tiers), fitter and products-only fit work; EuDetCode is external. **Missing:** per-step conservation monitor and error norms (WP 0.1), BKW/exp-6 tabulation (WP 0.3; a `BKW` EOS class exists in `src/thermodynamicModels` but is not tried). No committed harness report covers all tiers. |
| 1 0D physics | Partial | WP 1.1 met in the solver (isochoric only; isobaric missing). WP 1.2–1.5 are numpy models; the stage document marks them "PASSED", which overstates the position. WP 1.3 in the solver is at stage 1 of 2 (constant-$C_p$). Gate 1 is not met: no experimental ignition-delay comparison. |
| 2 1D numerics | Mostly complete | Gate 2 met (gunIB1D within 1.2%). WP 2.1 is weaker than reported: the CFD shock tube runs one Mach number with **non-dissociating** JANAF air (T2 1890 K) and is compared with the exact jump at that Mach number; the Mach 3–10 dissociation table is purely analytic. WP 2.4 uses Comp-B without Al (1720 kg/m³). Compressible drag not implemented. |
| 3 1D spherical | In progress | 3-phase charge, interface and Al cases run. Gate 3 is not met and is measured on the wrong case (C4 in blastFoam, 300 cells, 10–20% overpressure error) instead of the inert Comp-B/Al baseline. **Not started:** WP 3.2 transfer sweep (transfer is off), long-time runs, WP 3.4 impulse gain, WP 3.5 BKW trial. |
| 4 Dispersed frameworks | Not started | Only the existing Eulerian (blastEulerFoam) and Lagrangian (blastParcelFoam) solvers; no QMOM, no hybrid handover. |
| 5 Multi-D | Not started | 2D `ChargeCompBAl` tutorial set up. The first run crashed at 116 µs (floating-point exception); a restart from 110 µs passed that point and is at 126 µs. No `blastRadFoam` (radiation models exist in `src/radiationModels`). |
| 6 Validation | Not started | No experimental dataset identified. |

---

## 3. Findings

### 3.1 Physics gaps on the critical path

1. **Products → gas transfer is off in all Comp-B/Al cases.** `pressurePhaseChange` is configured in `fvModels.disabled` in `_reacting` and `_alCombustion` (since `124926e3` / `6f693f7b`) and in the three tutorials. This is the mechanism of report Section 3 and the precondition for phases (2) and (4). Neither the report nor `UPDATES.md` (Open) records that it is off.
2. **The "reacting" spherical case does not react.** It differs from the inert case only in `solveDragODE on` and `alphaCutoff 1e-3` on drag and heat transfer, plus laminar chemistry on a gas that contains only O2/N2. The verifier's chemistry check is that the file `2e-05/Qdot.gas` exists. This also explains the open "inert vs reacting U_Al edge difference" (2303 vs 1890 m/s): without the cutoff, drag accelerates the dilute leading edge ($\alpha \approx 10^{-3}$) of the inert case's Al.
3. **Al particle physics is not in the solver.** The CFD burn law is $r_b = a\,p^n$ (`src/diameterModels/reactionRates/pressureBased`). There is no ignition temperature, no melting, no oxidiser dependence beyond the `YoxidiserRef` scaling, and no Mach/Knudsen-dependent drag (none of the drag models in `interfacialModels/dragModels` depends on Mach number).
4. **High-temperature air.** Constant-$C_p$ extrapolation (stage 1) prevents the 5000 K cap but gives 14 400 K (spherical) and 24 800 K (planar). NASA-9 with three ranges and O2/N2/NO dissociation are needed (WP 1.3 stage 2).
5. **Time and length scales.** All Comp-B/Al cases stop at 40 µs in a 0.2 m domain (Z ≈ 0.2 m/kg^1/3 for the 0.95 kg charge). Kingery–Bulmash comparison, positive impulse and afterburning need 2–5 m and 1–5 ms (WP 3.1 specifies R = 5 m). No case has been run there.

### 3.2 Numerical robustness

- 2D `ChargeCompBAl`: floating-point exception in the products energy solve at 116 µs. It depends on the exact sequence of time steps (a restart passed it), so some cell sits on the edge of failing. Not analysed.
- `tutorials/blastEulerFoam/reactingParticles`: the time step collapses near 0.27 ms (not analysed).
- In the spherical cases the trace-gas temperature reaches the 50 000 K decode cap at initiation (9–45 of ~2540 steps).
- Verifiers do not check that the solver reached `End` (the runner does; scratch runs did not).

### 3.3 Documents that disagree with the code

| Location | Claim | Actual |
|---|---|---|
| `stage1_0d_unit_physics_kinetics.md` (note box and table) | WP 1.2–1.5 "Implemented & Verified", "PASSED", target `src/combustionModels/...`, `src/lagrangian/...` | numpy models only; the paths do not exist |
| `stage2_1d_canonical_multiphase.md` Section 3 | Houim 2 "0.02% error"; WP 2.1 "real-gas jump matched across Ma 3–10" | 0.02% was gunIB1D vs reference (blastEulerFoam vs gunIB1D: 1.20%); the CFD runs one Mach number without dissociation |
| `development_plan.md` Stage 0 row | "the case fits have not been regenerated" | regenerated 2026-10-13 |
| `development_plan.md` Stage 3 row | "`_interface` fails at 5 µs" | fixed 2026-10-14 (`alphaTrace`) |
| `development_plan.md` table headers | "Status (2026-10-10)", "as of `c3294dc1`" | rows updated to 2026-10-14 |
| `TBX-blastFoam.tex` Table `tab:regression_summary` | Tier 2 row `ChargeCompBAl_1D_AlInert`; "16 test cases" | tutorial removed (runner now uses blastEulerFoam `ChargeCompBAl_1D`); the runner has 18 tests |
| `TBX-blastFoam.tex` Summary, Stage 1 | "0.59% maximum pressure variation" | 3e-7 since 2026-10-13 |
| `TBX-blastFoam.tex` Stage 3 / table | `_reacting`: "finite-rate CH4/CO, $\dot Q_\text{gas}$ present" | no fuel in the gas (3.1, item 2) |

---

## 4. Changes to the original plan, and its weaknesses

### 4.1 Changes already made (justified by the work so far)

- **blastEulerFoam is the main solver** for aluminised charges. Under the original plan, blastFoam's single-velocity mixture would have been enough until Stage 4; the separate Al phase (report Section 2.1) needs per-phase momentum and energy.
- **Products-only eJWLC fit** with a separate inert Al phase ($\alpha_\text{Al} = 0.182$), instead of a full-mixture fit.
- **Conservative multi-fluid decode** and dilute-phase pressure/temperature equilibration (`alphaDilute`, stiff products/gas heat transfer), plus HLLC trace-phase impedance (`alphaTrace`). These were not in the plan and were needed to make phase (1) stable and conservative (`dilute_phase_conservation_plan.md`).
- **WP 1.3 extended** from 6000 K to 20 000 K (NASA-9, three ranges, dissociation in the solver).
- **Exit tolerances relaxed in practice.** For example, WP 1.1 asked for $|\Delta E|/E < 10^{-14}$, but the iterative decode gives ~5e-9. The plan should state attainable criteria.

### 4.2 Weaknesses of the original plan

1. **Counting "verified" without the solver.** Stage 1 accepts standalone numpy models as verification of physics the CFD does not contain. A WP should count as done only when the model runs in the solver and has a solver test.
2. **The Stage 6 decision matrix is pre-filled.** It already scores Hybrid at 92%, before any Stage 4 work. The scores should come from the Stage 4 benchmark.
3. **Stage 4 is too wide.** Building QMOM, dense DPM and a hybrid handover before there is a working afterburning baseline is costly. Selection should start with the two solvers that exist (Eulerian, Lagrangian) and add others only if a measured shortcoming needs them.
4. **No validation data identified.** Gate 5 needs experimental gauge data for a Comp-B/Al charge; the WP 1.4 ignition delay and WP 3.4 impulse gain (25–40%) also need sources. None are named in `report/external_sources` or the plans.
5. **Gate 3 reference.** Kingery–Bulmash is TNT data. The inert Comp-B/Al baseline needs a TNT-equivalence argument, or an aluminium-free Comp-B baseline, and a domain of ≥ 2 m.
6. **Sequential stages.** Phase transfer calibration (Stage 3) and Al ignition/combustion (Stage 1 solver work) do not depend on Stage 2, but the linear stage ordering hides that the critical path runs through them.
7. **Compute cost not planned.** At ~4 s per step for 160 000 cells (2D, 6 cores), multi-millisecond 2D/3D runs need AMR, more cores or larger steps. AMR is in Stage 5 but is needed earlier for the long 1D/2D runs.

---

## 5. Step-by-step plan

Ordered by the critical path to Gate 5. Effort: S ≈ days, M ≈ 1–3 weeks, L ≈ > 3 weeks of focused work (indicative). Each step lists its exit criterion.

### Phase A: Consolidate (S–M)

| # | Step | Exit criterion |
|---|---|---|
| A1 | Correct the documents in Section 3.3; mark numpy-only WPs as "model only"; record in `UPDATES.md` that phase transfer is off in the Comp-B/Al cases | Docs match the code |
| A2 | Verifiers check `End` in the solver log; the `_reacting` verifier checks heat release and product species, not file existence | A crashed or non-reacting run fails |
| A3 | Run the full harness (`--tier all`) and commit `blast_test_report.md` | One report covers all 18 tests |
| A4 | Finish the 2D `ChargeCompBAl` run; find the failing cell of the 116 µs crash (restart from the 110 µs write with the original step sequence) and fix it | 2D tutorial runs to 200 µs reproducibly |
| A5 | Remove the 50 000 K trace-gas cap hits at initiation (initial trace-gas state in the charge) | No cap hits |

Progress (2026-10-14): A1 done (stage documents, development plan, report, `UPDATES.md`). A2 done: the verifiers check `End`, and the `_reacting` verifier checks for heat release, so it fails until B1. A5 partly done: one source (HLLC flux from a residual upwind side with `alphaTrace`) fixed, cap hits 11 → 3 (inert) but 45 remain in `_alCombustion` (`UPDATES.md`, Open). A3 done: full run of all 18 tests on the final code, 15 pass; `_reacting` and the (now Mach 14.5) air shock fail by design, gunIB1D passes when not cut off by the runner timeout (raised). A4: the 2D crash was a sign error in the upstream fluid-phase residual stabilisation; fixed, the run reaches 200 µs from a 170 µs restart; a full run from 0 is still to be made.

### Phase B: Products → reactive gas in the spherical cases (M) — critical

| # | Step | Exit criterion |
|---|---|---|
| B1 | Enable `pressurePhaseChange` (products → gas, freeze-out composition) in `_reacting` and `_alCombustion`; check global mass/energy | Conservation to ~1e-8; no pressure spike at transfer |
| B2 | Extend the domain to R = 2 m, end time ≥ 2 ms (AMR or graded mesh) | Runs complete; shock at 2 m captured |
| B3 | WP 3.2 sweep: $p_\text{threshold} \in \{50, 100, 200, 500\}$ bar, $\tau \in \{10^{-7} \ldots 10^{-4}\}$ s | Gauge $p(t)$ variance < 1.5% for 100–200 bar, or a documented sensitivity |
| B4 | Product afterburning: CO/H2 + air in the reacting case (current 1-step reactions; then a reduced CO/H2/O2 mechanism) | $\dot Q_\text{gas}$ in the mixing layer; afterburn energy vs CO/H2 inventory |

Progress (2026-10-14): B1 done for the spherical `_reacting` and `_alCombustion` cases (transfer below 100 bar, τ = 20 µs, `seulex` chemistry); the products burn with the air and the `_reacting` test passes. The transfer adds no visible energy error, but a closed-domain check found that total energy is not conserved with a particle phase (+2.8% spherical, −0.55% planar gas/particle shock tube), independent of the transfer. This has to be fixed before the B3 sensitivity study and E2 impulse gain can be trusted; added as A6 below.

| A6 | Find and fix the energy error of the particle/fluid coupling (pressure-gradient work, granular energy exchange); closed-domain energy check in the test suite | Total energy conserved to ≤ 1e-6 in `gunib1d_houim2` and the spherical cases while closed |

Progress (2026-10-14): A6 done for the coupling: three energy errors fixed (drag ODE without energy, missing pressure work on the particle-occupied volume, implicit stage rebuilding the conserved fields from bounded primitives). Houim 2 conserves energy to 6e-7 (target met), the spherical cases to 4–7e-4 (remaining: end of the programmed burn). Energy checks are in the verifiers. All 18 tests pass except the air shock (fails by design).

Progress (2026-10-15): B2 done: the 2 m / 2 ms run of the `_alCombustion` setup is the tutorial `tutorials/blastEulerFoam/ChargeCompBAl_1D_spherical_2m` (graded mesh, 1 h 45 min serial, too long for the test runner). B3 and B4 are open.

### Phase C: High-temperature air (M) — WP 1.3 stage 2

Extended 2026-10-14 to local chemical equilibrium (the one-step finite-rate chemistry fails at high temperature); detailed plan in `phaseC_highT_equilibrium_plan.md`. Done 2026-10-14 as steps C1–C7 of that plan (the table below is the original outline).

| # | Step | Exit criterion |
|---|---|---|
| C1 | NASA-9 three-range thermo class (200–1000–6000–20 000 K), instantiated for the multicomponent gas; species N2, O2, NO, N, O, Ar | $c_p$, $h$, $s$ match CEA to 0.2% to 20 000 K (0D solver test) |
| C2 | Dissociation reactions O2 ⇌ 2O, N2 ⇌ 2N, NO with third bodies (Park rates); stiff chemistry check | 0D equilibrium composition matches CEA |
| C3 | WP 2.1 in the solver: shock tube at Ma 3, 5, 8, 10 with dissociating air | Jump conditions within 2.5% (CFD) of equilibrium Rankine–Hugoniot |
| C4 | Re-run the spherical and planar cases | Shocked-air T ≤ ~10 000 K; update numbers |

### Phase D: Al particle physics in the solver (L) — WP 1.2, 1.4, 1.5, 2.3

| # | Step | Exit criterion |
|---|---|---|
| D1 | Al thermo with melting latent heat (933 K, 397 kJ/kg) in the granular phase | 0D heating test shows the plateau |
| D2 | Mach/Knudsen-dependent drag (Henderson or Loth 2008/2021) and a Nusselt correlation with compressibility | 0D/1D single-particle tests vs correlations; gunIB1D Houim cases unchanged |
| D3 | Ignition criterion ($T_p \ge T_\text{ign}$) and Beckstead $d^n$ multi-oxidiser rate as a new `diameterReactionRate` | 0D bomb: burn time vs Beckstead correlation within 5% for 1–50 µm |
| D4 | Al2O3 condensation (gas → condensed oxide) with latent heat; check the Al2O3(s) species route already in the gas | Stoichiometric Al/air 0D flame locks near 3800 K; total energy 31 MJ/kg Al |
| D5 | Replace the numpy verifiers by solver 0D tests (keep numpy as the reference) | Tier 1 tests use blastEulerFoam |
| D6 | Collect experimental data: ignition delays, burn times (Beckstead, Gurevich, shock-tube data) | Gate 1 comparison |

Progress (2026-10-15): D1 done (`eMelting` thermo; `tests/0D/alMelting_0D`), D2 done (`Loth` drag, `tests/0D/particleDrag_0D`; `Kavanau` Nusselt number, `tests/0D/particleHeating_0D`; the pair Prandtl number used c_v, fixed); D3 done (`Beckstead` diameter rate with `Tign`; `tests/0D/alBurning_0D`). D3 found two errors in `reactingDiameterModel`: the number density used the initial diameter (a burnt-out cloud kept 1/e of its mass) and a volume floor of 1e-15 m³ (micron particles lost mass up to 2000 times too slowly); fixed, and the reset volume fraction is now the entry `alphaMin` (default 1e-4, which stops burning in the dilute far field). D4 (2026-10-15): Al2O, AlO2, Al2O2 and liquid Al2O3 in the gas species data; the equilibrium of stoichiometric Al(s)/air at constant volume is 4070 K at 1 atm air density (12.9 bar) and 4561 K at 10 times, with most of the Al as liquid Al2O3, against an independent solution (`tests/0D/verify_al_flame.py`); 31.05 MJ/kg Al. Open: D5 (the numpy reference models stay), D6, and using the new models and species in the Comp-B/Al cases (E2).

### Phase E: 1D spherical validation (M)

| # | Step | Exit criterion |
|---|---|---|
| E1 | Aluminium-free Comp-B spherical baseline, R ≥ 2 m, vs Kingery–Bulmash with TNT equivalence | Gate 3: peak overpressure within 5% for Z ≥ 1 |
| E2 | Inert vs reacting Comp-B/Al with B, C, D in place | WP 3.4: impulse gain quantified; Al burnout radius/time |
| E3 | Optional WP 3.5: BKW single-EOS trial (existing `BKW` EOS) on the products | Decision note: keep eJWLC → gas or not |

Progress (2026-10-15): E1 done for TNT and aluminium-free Comp-B (`tests/1D_spherical/tnt_blast`, `compB_blast`, `verify_blast_kb.py`, reference `kingery_bulmash.py`): Gate 3 met for Z ≥ 4 at 1 mm (TNT within 3.4%) but not at Z 1.5–3 (5–8% low); the overpressure depends on the resolution of the charge (4 mm: 5–15% low), not on the limiter; the positive impulse is 14–17% low for TNT and 5–9% for Comp-B independent of the mesh, consistent with the missing afterburning. B3 (preliminary, 0.2 m / 40 µs `_alCombustion`, 16 runs): the gauge peaks at 80 / 100 / 150 mm vary by 0.07 / 0.56 / 2.9% for transfer pressures 100–200 bar and τ 0.1–100 µs; three runs with high transfer pressure and fast transfer (≥ 200 bar with τ 0.1 µs, 500 bar with τ ≤ 1 µs) aborted in the gas temperature inversion: the transfer acted in the burn front, where λ > 0 before the pressure has risen, and moved unreacted charge into the gas. Fixed by `lambdaMin` (default 0.99) in `pressurePhaseChange`; all 16 runs complete. For 100–200 bar and τ ≥ 1 µs the peaks vary by 0.05 / 0.56 / 2.9% at 80 / 100 / 150 mm; a near-instant transfer (τ = 0.1 µs at 200 bar) lowers the 100 mm peak by 5%. The B3 criterion should be evaluated on the 2 m case.

Progress (2026-10-19): E2 on the 2 m / 2 ms case (`tutorials/blastEulerFoam/ChargeCompBAl_1D_spherical_2m/setup_e2.py`, `compare_e2.py`; Phase D Al models, burning to oxide with the gas-phase oxidisers, transported ignited fraction, 28-species gas; the Al couples to the products by drag and heat since the blending fix of 2026-10-16). Products afterburning: +2–5% peak, +3–7% impulse beyond 0.3 m. Al: +10% peak at 0.5 m, +36–52% from 1 m, impulse at 1 m doubled; 90–95% burns with air O2 at 1–1.5 m ahead of the products. Ignition at 2050 K (intact shell), 932.5 K (shell ruptured by the detonation) and 0: 66.5 / 80.4 / 80.4% burnt at 1.9 ms, peaks within 1%; the products melt all the Al within 20 µs, so the ruptured-shell case ignites at once. In 1D the Al contribution is probably overestimated (no mixing limit). WP 3.4 impulse gain quantified for 1D; burnout radius not reached by 2 ms (80% burnt). B3 on the 2 m gauges (products only): 200 bar within 0.3%; a fast transfer (2 µs) changes the peaks by up to 2.3% and the impulse at 0.3 m by 10%, so the 1.5% criterion holds for 100–200 bar at τ = 20 µs.

### Phase F: Dispersed-phase method (M–L) — narrowed Stage 4

| # | Step | Exit criterion |
|---|---|---|
| F1 | Same 1D/2D Al cloud problem in blastEulerFoam (Eulerian) and blastParcelFoam (Lagrangian, with the Phase D models) | Dispersion, burn rate, cost compared on identical cases |
| F2 | Only if F1 shows a clear deficit (e.g. polydispersity, edge diffusion): QMOM or hybrid handover | Decision matrix scored from F1 data (Gate 4) |

### Phase G: Multi-D (L) — Stage 5

| # | Step | Exit criterion |
|---|---|---|
| G1 | WP 5.1: 2D axisymmetric repeat of the 1D spherical case with AMR | Angular variance < 0.5%; overlay on 1D |
| G2 | WP 5.5: AMR sensors and cost (needed for G3/H1) | > 75% saving vs uniform fine mesh |
| G3 | Cylindrical charges (L/D) and mixing (RMI/RTI, LES vs RANS, PaSR) | Mixing-layer growth vs literature; afterburn sensitivity to turbulence closure |
| G4 | Radiation (existing `radiationModels`, P1) where alumina smoke matters | Radiative loss fraction quantified |

Progress (2026-10-20): G1 accepted with revised criteria (`tests/2D_axi/ChargeCompBAl_2D_axi_spherical`, AMR from 2 mm to 0.5 mm, 6 processes). The mean of the 30–60° rays overlays the 1D case within 0.74% in arrival and 2.7% in peak; between 10 and 80° the peaks vary by up to 3% and the arrival by 1.5%, against a 1–3% change of the 1D case itself at 0.25 mm. The < 0.5% target is not met: the HLLC flux has a carbuncle at the axis and the symmetry plane (peaks up to 9% low, a trace of products running along both to the shock), also on the uniform mesh. A hybrid HLLC/HLL flux (`hybridHLL`, pressure-ratio sensor) removes it (axis and plane within the range of the other directions) but raises the air-shock peaks by 1–3%; made the default on 2026-10-20 (full suite passes; overlay tolerances of the 2D test set to the 1D resolution sensitivity, arrival 1.5%, peak 5%). AMR: 2.2 h against 5.7 h on the uniform mesh (62% saving; G2 target 75%). Also fixed: `temporalMax`/`temporalIntegral` maps were wrong after load balancing.

Progress (2026-10-21): G2 met (`tutorials/blastEulerFoam/ChargeCompBAl_2D_axi_2m`, the full E2 setup on a 2 m wedge, 6 processes). The run starts from the 1D solution mapped at 60 µs (`rotateFields -refine` in parallel after `decomposePar`, which does not carry the refinement history), refines on 20% jumps of density and pressure and at the Al cloud edge (`scaledDelta`; `densityGradient` kept the whole fireball refined), and lowers the finest level with time like the 1D graded mesh (2 / 4 / 8 mm). To 2 ms: 5.6 h, 17 700–24 900 cells, 3.7% of the cell-steps of a uniform mesh at the same levels (96% saving). Against 1D at 0.3–1.9 m: arrival within 1.1%, peak −0.7 to +9%, positive impulse +1 to +27% (0.5 m); the products afterburn more in 2D (14% less CO and 12% less H2 at 2 ms), the Al burns as in 1D (60 / 57%). The extra afterburning is numerical mixing: the refinement criteria lose the products/air mixing layer after 0.2 ms (from 0.5 ms all its cells are at 16 mm; layer 3–4 times wider than in 1D), so a criterion that follows the mixing layer is needed before G3. Angular spread of the peaks 1–5%, axis and symmetry plane within the other directions. Fixed on the way: cells above a decreasing `maxRefinement` were not coarsened (`polyMeshHexRefiner`), the multi-estimator maximum level (`errorEstimator`), and the detonation energy counted twice on restarts (`activationModel::initESource`; it affected the 170 µs restart of the 2D tutorial).

### Phase H: Validation (L) — Stage 6

| # | Step | Exit criterion |
|---|---|---|
| H0 | **Now:** identify the experimental Comp-B/Al dataset (charge mass, geometry, gauge distances) | Dataset named in the plan |
| H1 | 3D or 2D-axi run of the trial geometry, inert and reacting | Gate 5: arrival ±3%, peak ±7%, impulse ±8% |
| H2 | Production guidelines and tutorials | WP 6.4 |

### Order and dependencies

```
A (consolidate) ──> B (transfer) ──> E1/E2 (1D validation) ──> G (multi-D) ──> H1
        │                  ^                    ^
        ├──> C (high-T air)┘                    │
        └──> D (Al physics) ────────────────────┘──> F (method choice)
H0 (find data) and D6 (kinetic data) start now
```

The near-term order: A1–A3, then B1 together with C1, then D1–D3. B and C are independent, and D can start in parallel once A is done.
