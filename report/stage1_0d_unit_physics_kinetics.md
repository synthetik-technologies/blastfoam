# Stage 1: 0D Unit Physics & Kinetics Verification

**Goal**: Isolate and rigorously verify individual chemical kinetic mechanisms, single-particle thermal ODEs, aluminium ignition thresholds, high-temperature gas dissociation thermodynamics, and alumina condensation latent heat in zero-dimensional configurations.

---

## 1. Overview & Work Packages

```mermaid
graph TD
    WP11["WP 1.1: 0D-1 Phase Transfer & Energy Conservation"] --> WP12["WP 1.2: 0D-2 Al Oxidation Multi-Oxidizer Kinetics"]
    WP11 --> WP13["WP 1.3: 0D-3 High-T Air Dissociation & Thermodynamics"]
    WP12 --> WP14["WP 1.4: 0D-4 Particle Thermal ODE, Melting & Ignition"]
    WP12 --> WP15["WP 1.5: 0D-5 Alumina Condensation & Latent Heat Coupling"]
```

> [!NOTE]
> **Status (2026-10-14).** Only WP 1.1 is implemented and verified in the solver. WP 1.2–1.5 and the dissociation part of WP 1.3 are standalone numpy models, not solver code; they are reference models for the solver implementation (see `status_review_2026-10-14.md`, Phase C and D).
> - **WP 1.1**: blastEulerFoam [`tests/0D/phaseTransfer_0D/`](file:///opt/OpenFOAM/blastfoam/tests/0D/phaseTransfer_0D) and `phaseTransfer_0D_3phase` (inert Al), verified by [`tests/0D/verify_0d_phase_transfer.py`](file:///opt/OpenFOAM/blastfoam/tests/0D/verify_0d_phase_transfer.py): total energy conserved to 4.2e-9 / 2.2e-9 (the iterative decode limits this; the $10^{-14}$ target below is not attainable), pressure constant to 3e-7, > 99% conversion. Isochoric only; the isobaric case is missing.
> - **WP 1.2, 1.4, 1.5**: numpy model [`tests/0D/verify_al_particle_kinetics.py`](file:///opt/OpenFOAM/blastfoam/tests/0D/verify_al_particle_kinetics.py) (melting plateau, ignition at 2050 K, Beckstead multi-oxidiser rate, alumina condensation). Not in the solver, which uses the empirical `pressureBased` burn rate; no comparison with experimental data.
> - **WP 1.3**: numpy model [`tests/0D/verify_highT_air_thermo.py`](file:///opt/OpenFOAM/blastfoam/tests/0D/verify_highT_air_thermo.py) to 6000 K. In the solver: JANAF (7-coefficient) to 5000 K, extrapolated with constant $c_p$ above (no dissociation) since 2026-10-14.
> - Tier 1 of [`bin/run_blast_tests.py`](file:///opt/OpenFOAM/blastfoam/bin/run_blast_tests.py) runs all of these.

| Work Package | Description | Target Code / Dictionaries | Verification Case | Status |
|---|---|---|---|---|
| **WP 1.1** | 0D-1 Phase transfer verification (Isobaric & Isochoric) | `src/thermodynamicModels/fvModels/pressurePhaseChange/` | `tests/0D/phaseTransfer_0D/`, `phaseTransfer_0D_3phase` | **Solver: isochoric passed**; isobaric missing |
| **WP 1.2** | 0D-2 Multi-oxidizer aluminium reaction kinetics ($d^n$) | solver: not implemented (`src/diameterModels/reactionRates/` has `pressureBased`, `Arrhenius`, `constant`, `surfaceReactionRate`) | `tests/0D/verify_al_particle_kinetics.py` | **Model only** (numpy) |
| **WP 1.3** | 0D-3 High-temperature air dissociation ($T \le 6000\,\text{K}$) | `src/thermodynamicModels/specie/thermoModels/janaf/` | `tests/0D/verify_highT_air_thermo.py` | **Model only** (numpy, 6000 K); solver stage 1 (constant $c_p$) done |
| **WP 1.4** | 0D-4 Particle heating, melting plateau, oxide shell ignition | solver: not implemented | `tests/0D/verify_al_particle_kinetics.py` | **Model only** (numpy) |
| **WP 1.5** | 0D-5 Post-combustion $\text{Al}_2\text{O}_3$ phase change & latent heat | solver: not implemented | `tests/0D/verify_al_particle_kinetics.py` | **Model only** (numpy) |

---

## 2. Work Package Details

### WP 1.1: 0D-1 Phase Transfer Isochoric & Isobaric Verification
- **Formulation**:
  - Test closed single-cell system with `pressurePhaseChange` fvModel transitioning from `compBAl` (`eJWLC`) to `gas` (`Noble-Abel`).
  - Mass transfer is triggered over a selected pressure interval $[p_\text{lower}, p_\text{upper}]$ (e.g. $[50, 150]\,\text{bar}$) rather than temperature intervals due to JWL temperature inaccuracies.
  - Verify constant volume ($V = \text{const}$):
    $$\frac{d(\rho e_{\text{total}})}{dt} = 0, \quad \rho_{\text{final}} = \rho_{\text{initial}}, \quad p_{\text{final}} \approx p_{\text{initial}}$$
  - Energy conservation:
    $$\Delta E / E_0 = 0 \text{ (machine precision)}$$
- **Configuration**:
  - Sweep transfer timescale $\tau \in [10^{-8}, 10^{-4}]\,\text{s}$ and pressure smoothing width $p_{\text{range}}$.
- **Exit Criteria**:
  - Exact total energy conservation ($|E(t) - E(0)| / E(0) < 10^{-14}$).
  - Maximum pressure deviation $\Delta p / p < 0.01\%$ across the entire transition.

### WP 1.2: 0D-2 Multi-Oxidizer Aluminium Kinetics & Beckstead Bomb Benchmark
- **Formulation**:
  - Heterogeneous reaction model with oxidizers $\text{O}_2, \text{H}_2\text{O}, \text{CO}_2$:
    $$4\text{Al} + 3\text{O}_2 \to 2\text{Al}_2\text{O}_3, \quad 2\text{Al} + 3\text{H}_2\text{O} \to \text{Al}_2\text{O}_3 + 3\text{H}_2, \quad 2\text{Al} + 3\text{CO}_2 \to \text{Al}_2\text{O}_3 + 3\text{CO}$$
  - Burn rate: $\frac{d(d_p)}{dt} = -\frac{2 K_{\text{comb}}}{\rho_p d_p^{n-1}} X_{\text{eff}}^{0.9}$ with $X_{\text{eff}} = X_{\text{O}_2} + 0.6 X_{\text{H}_2\text{O}} + 0.22 X_{\text{CO}_2}$.
- **Exit Criteria**:
  - Burn time $t_b(d_p)$ matches Beckstead experimental correlations for $d_p \in [1\,\mu\text{m}, 50\,\mu\text{m}]$ within 5% across varied $\text{O}_2/\text{H}_2\text{O}/\text{CO}_2$ mixtures.

### WP 1.3: 0D-3 High-Temperature Air Dissociation & Thermodynamics
- **Formulation**:
  - Verify NASA-9 polynomials for air species ($\text{N}_2, \text{O}_2, \text{NO}, \text{N}, \text{O}, \text{Ar}$) up to $6000\,\text{K}$.
  - Validate forward and backward dissociation reaction rates: $\text{O}_2 + M \rightleftharpoons 2\text{O} + M$, $\text{N}_2 + M \rightleftharpoons 2\text{N} + M$, $\text{NO} + M \rightleftharpoons \text{N} + \text{O} + M$.
- **Exit Criteria**:
  - Equilibrium species mole fractions and specific heats ($c_p, c_v, \gamma$) match CEA/JANAF database within 0.2% for $T \in [300, 6000]\,\text{K}$.
- **Status / extension (2026-10-14)**: the solver gas phases use 7-coefficient JANAF fits to 5000 K, since 2026-10-14 extrapolated with constant $c_p$ above that (no dissociation). Shocked air in the Comp-B/Al cases reaches 14 400 K (spherical) and 24 800 K (planar), so the target range is extended to 20 000 K: NASA-9 polynomials with three temperature ranges (200–1000, 1000–6000, 6000–20 000 K) in a C++ thermo, and the O2/N2/NO dissociation reactions above in the gas-phase chemistry, verified against CEA to 20 000 K.

### WP 1.4: 0D-4 Particle Thermal ODE, Melting Plateau & Ignition Delay
- **Formulation**:
  - Single particle energy equation:
    $$\rho_p c_p \frac{V_p}{A_p} \frac{dT_p}{dt} = h_{\text{conv}}(T_g - T_p) + \epsilon \sigma (T_{\text{rad}}^4 - T_p^4) - \dot{q}_{\text{loss}}$$
  - Latent heat of melting: plateau at $T_{\text{melt}} = 933\,\text{K}$ until $\int \dot{q} \, dt = m_p \Delta H_{\text{fus}}$ ($\Delta H_{\text{fus}} = 397\,\text{kJ/kg}$).
  - Ignition activation when $T_p \ge T_{\text{ign}} \approx 2050\text{--}2300\,\text{K}$ (oxide shell melting/rupture).
- **Exit Criteria**:
  - Thermal response reproduces experimental laser-heated particle ignition delay curves across heating rates from $10^5$ to $10^8\,\text{K/s}$.

### WP 1.5: 0D-5 Post-Combustion $\text{Al}_2\text{O}_3$ Condensation & Enthalpy Buffering
- **Formulation**:
  - Enthalpy evolution coupling condensation $\Delta H_{\text{cond}} \approx 18.6\,\text{MJ/kg}$ and solidification $\Delta H_{\text{fus}} \approx 1.1\,\text{MJ/kg}$:
    $$\text{Al}_2\text{O}_{3(\text{g})} \to \text{Al}_2\text{O}_{3(\text{l})} \to \text{Al}_2\text{O}_{3(\text{s})}$$
- **Exit Criteria**:
  - Adiabatic flame temperature of stoichiometric $\text{Al} + \text{Air}$ mixture locks at alumina condensation plateau ($\approx 3800\,\text{K}$ at blast pressures) and releases exactly $31.0\,\text{MJ/kg}$ Al total combustion energy.
