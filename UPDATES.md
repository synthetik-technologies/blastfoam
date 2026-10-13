# Library Updates (Compared to `origin/develop` @ 71b97fd)

This document summarises the updates, new features, thermodynamic models, numerical changes and bug fixes in this fork relative to `origin/develop` (commit `71b97fd7c9195fdb342ed4230d258dd1dde80b55`, the last merge of upstream `synthetik-technologies/blastfoam`).

*Last synchronised with the code on 2026-10-10 (after `9efec4f1`, including the fixes in Section 7).* Commits covered: `4cb49946`, `0a181853`, `124926e3`, `6f693f7b`, `4fd73e06`, `2971df46`, `68cbb94b`, `bdccfc0e`, `8b274c09`, `fa2c2283`, `c3294dc1`. See [Known issues](#7-known-issues) for open defects found while reviewing these changes.

---

## 1. New Thermodynamic Models & Equations of State

### 1.1 Multicomponent & Reacting Phase Fluid Thermodynamics
- **`multicomponentPhaseFluidBlastThermo`** (`src/thermodynamicModels/derived/multicomponentFluidBlastThermo/`):
  - Multicomponent phase fluid thermodynamics for multiphase systems (`phaseFluidBlastThermo`), so individual phases can consist of multiple chemical species with species transport and mixture properties.
- **`reactingPhaseFluidBlastThermo`** (`src/thermodynamicModels/derived/reactingPhaseFluidBlastThermo/`):
  - Pairs phase multicomponent thermodynamics with `blastCombustionModel` for per-phase chemical reactions.
  - Writes the phase heat release rate `Qdot.<phase>` ($[\text{W/m}^3]$, `AUTO_WRITE`), set from `combustion_->Qdot()` in `correct()`.

### 1.2 `eJWLC` Equation of State
- **Energy-based Mie-Grüneisen JWL with 3-term isentrope** (`src/thermodynamicModels/specie/equationOfStates/equationOfStates/eJWLC/`):
  - Full 3-term ($A, B, C, R_1, R_2, \omega$) isentrope pressure formulation:
    $$p(\rho, e) = A \exp(-R_1 V) + B \exp(-R_2 V) + C V^{-(\omega+1)} + \omega \rho \left(e - e_s(V)\right)$$
    where $V = \rho_0 / \rho$ and $e_s(V)$ is the isentrope energy (plus an optional constant `E0`).
  - Analytic CJ-isentrope reference temperature, with $V_\text{CJ} = \rho_0/\rho_\text{CJ}$:
    $$T_s(V) = T_\text{CJ} \left(\frac{V_\text{CJ}}{V}\right)^\omega$$
  - Requires `Cv`, `T_CJ` and `rho_CJ` in the dictionary. Gives isentropic and off-isentropic temperatures when coupled with `eConst` caloric models.
  - Optional `omegaT` (default `omega`, added 2026-10-10) is the exponent of $T_s(V)$. It changes only the energy–temperature relation, so $T_s$ can pass through both $T_\text{CJ}$ and the gas-EOS temperature at the transfer state without changing the pressure: $\omega_T = \ln(T_\text{CJ}/T_\text{ref})/\ln(V_\text{ref}/V_\text{CJ})$.
  - Registered in `src/thermodynamicModels/specie/include/forBlastGases.H` and `forDetBlastGases.H`.

### 1.3 Solid Equations of State
- **`forBlastSolids.H`**: Added `Murnaghan`, `BirchMurnaghan2` and `BirchMurnaghan3` for solid thermodynamics.
- **`rhoConst`**: Added `rho(p, T)` and `rho0()` accessors.
- **`basicSolidBlastThermo::updateRho(p)`**: New; used by the compressible granular phase (Section 3.3).

---

## 2. Finite Volume Models (`fvModels`) & Phase Change

### 2.1 Pressure-Driven Phase Change Model (`pressurePhaseChange`)
Location: `src/thermodynamicModels/fvModels/pressurePhaseChange/`, built into `libblastThermodynamics` (moved from `src/finiteVolume` in `c3294dc1`, see Section 5). Also selectable under the alias `pressurePhaseTransfer`.

- Pressure-induced mass transfer from a source phase to one or more multicomponent target phases, driven by a pressure interval instead of a temperature threshold. This avoids the known inaccuracy of empirical JWL temperature/caloric formulations.
- Interval input (first match wins): `pInterval (pMin pMax)`, `pMin`/`pMax`, `pLower`/`pUpper`, or `pThreshold` (upper bound; `pLower` optional, default 0). A legacy `Tthreshold` entry gives a warning and falls back to $[0, 100\,\text{bar}]$.
- Transition factor, for $p_\text{lower} > 0$:
  $$\Gamma(p) = \tfrac{1}{4}\left(1 + \tanh\tfrac{p_\text{upper} - p}{p_\text{range}}\right)\left(1 + \tanh\tfrac{p - p_\text{lower}}{p_\text{range}}\right)\,H\!\left(p_\text{upper} + 2p_\text{range} - p\right)\,H\!\left(p - p_\text{lower} + 2p_\text{range}\right)$$
  where $H$ is a hard step (`pos`). For $p_\text{lower} \le 0$ the one-sided form $\tfrac{1}{2}\left(1 + \tanh\tfrac{p_\text{upper} - p}{p_\text{range}}\right)H(p_\text{upper} + 2p_\text{range} - p)$ is used.
  - $p_\text{range}$ defaults to $\max(0.05\,(p_\text{upper} - p_\text{lower}),\,1\,\text{Pa})$; an explicit `pRange` (or `Prange`) overrides it.
- Transfer rate $\dot m = \Gamma\,\alpha_s \rho_s/\tau$ (`tau`, default $10^{-4}$ s). Since `0a181853`, if a `lambda.<sourcePhase>` field exists (programmed burn), $\dot m$ is multiplied by it in the mass-source path.
- `targetFractions` sets the target-phase species distribution, as mass fractions by default. With `isMoleFraction true` they are converted on first use to mass fractions, $Y_i = X_i W_i/\sum_j X_j W_j$, using the molecular weights of the target-phase thermo. The converted values are printed to the log. A target specie that is not in the phase, or a target phase without multicomponent thermo, is a fatal error.
- Energy and momentum: all phases share an energy datum, reconciled through the products `E0` (offline, `utilities/thermo_fitter`). An optional constant datum shift `eOffset` [J/kg] (default 0) can be added to the transferred mass. No state-dependent offset is computed (see Section 7).
  - blastEulerFoam (multi-fluid, per-phase `alphaRho`, `alphaRhoU`, `alphaRhoE`; since 2026-10-10): the transferred mass carries the source-phase velocity $\mathbf{u}_s$ and specific total energy $E_s = (\alpha\rho E)_s/(\alpha\rho)_s$. The source loses $\dot m(1, \mathbf{u}_s, E_s)$; each target gains $f_k\dot m(1, \mathbf{u}_s, E_s + e_\text{off})$. Volume moves with the mass ($\pm\dot m/\rho_s$), so internal (not enthalpy) energy is carried. Mass, momentum and datum-corrected energy are conserved exactly.
  - blastFoam (single mixture energy and velocity): `rhoE` source $\dot m\,e_\text{off}$ only.
- Note: the code was introduced in `4cb49946` as a new model. No `temperaturePhaseChange` model was ever committed to this repository's history.

---

## 3. Solver & Compressible System Changes

### 3.1 `fvModels` Integration in Compressible Systems
- **`compressibleBlastSystem`**: `addSources()` adds energy source terms from active `fvModels` (`models().source(rhoEff(), e_)`).
- **`multiphaseCompressibleSystem`**: `solveMass()` adds phase source terms from active `fvModels` to both `alpha` (`models().source(alphas_[i])`) and `alphaRho` (`models().source(alphas_[i], rhos_[i])`).

### 3.2 Two-Phase System Stability & Reconstruction
- **`twoPhaseCompressibleSystem::decode()`**:
  - Volume fraction cutoff ($10^{-15}$) to prevent EOS evaluation in empty cells.
  - For phase volume fractions $\alpha < 0.01$, phase density is capped at a hard-coded $5000\,\text{kg/m}^3$.
  - The `transportPhaseDensity_` branch and commented-out phase-change code were removed.

### 3.3 `blastEuler` Module / `blastEulerFoam` (multiphase Euler–Euler)
Files under `applications/modules/blastEuler/` and `applications/solvers/blastEulerFoam/`.
- **Courant number** (`4cb49946`): phase Courant numbers are $\alpha$-weighted, and the reported `CoNum` is the summed $\alpha$-weighted value instead of the maximum phase Courant number (`blastEulerSolver.C`, `regionSolver/blastEulerRegionSolver.C`, `EigenCourantNos.H`).
- **Mass-exchange Courant number** (`8b274c09`): uses only net mass loss, $\max(-\sum\dot m, 0)\,\Delta t/(\alpha\rho)$, reduced with `gMax` (`EigenCourantNos.H`).
- **Drag** (`interfacialModels/dragModels/dispersedDragModel`, `dragODE`):
  - (`4cb49946`) $K_i$ scales as $d^2/\max(d, d_0)^4$, with $\nu \ge 10^{-15}$.
  - (`6f693f7b`) New `alphaCutoff` (falls back to `alphaMin`, default 0): drag is zeroed where either phase is below $\max(\alpha_\text{residual}, \alpha_\text{cutoff})$. `dragODE` skips cells with no active drag.
- **Heat transfer** (`RanzMarshall`, `spherical`): $K$ scales as $d/\max(d, d_0)^3$ instead of $1/d^2$. (2026-10-10) Fixed a dangling reference to the `tmp` returned by `d()`.
- **Kinetic theory**:
  - `kineticTheorySystem` collisional `lambda` term rewritten per cell with residual-$\alpha$ guards.
  - Bounded denominators in the `CarnahanStarling` and `SinclairJackson` radial models.
  - $\sqrt{\max(\Theta, \text{small})}$ in Gidaspow/Syamlal conductivity and viscosity.
  - **Syamlal viscosity: the collisional term expression was changed** (a formula change, not only a guard).
- **`phaseSystem`**:
  - Velocity, pressure–temperature-equilibrium and temperature relaxation exchange are zeroed where either phase is residual.
  - Mixture temperature averages only active phases (fallback 288.15 K); mixture density is limited.
  - `decode()` decodes granular phases first, rescales the carrier volume fraction, then decodes fluid phases.
  - (2026-10-12) `equilibrateDiluteVolumes()`, called in `decode()`: dilute fluid phases ($\alpha <$ `alphaDilute`, default 0.01, 0 disables) are brought towards pressure equilibrium with the carrier fluid phase by exchanging volume. It uses a damped Newton iteration in which the $\bar p\,\Delta\alpha$ work of each increment moves between the two phase energies. Only volume fractions and total energies change, so mass, momentum and energy are conserved. The dilute phase is compressed by at most a factor of two per call, and the energies stay within the thermo temperature ranges. See `report/dilute_phase_conservation_plan.md`.
- **`fluidPhaseModel` decode**:
  - (`4cb49946`) Density blending towards `residualRho` below $\alpha = 10^{-3}$, velocity blending to zero below $10^{-2}$, and an energy/temperature reset to 288.15 K below $\alpha = 0.01$, all written back to the conserved variables (non-conservative). `solveImplicit` calls `decode()` instead of `thermo.correct()`.
  - (2026-10-12) Made conservative: `alphaRho`, `alphaRhoU` and `alphaRhoE` are only modified where the phase is residual ($\alpha \le$ `residualAlpha`, reference state 288.15 K), as upstream. The velocity cap (20000 m/s) and the energy clamp to the thermo temperature range act on the primitive variables only.
  - (2026-10-13) The decoded density is no longer bounded below by `residualRho` (upstream `rho.max(residualRho)`); `residualRho` is only the reference density of residual cells. The gas `residualRho 1.2` held expanded gas at ambient density in up to 90 cells of the spherical runs and overstated its pressure up to twentyfold.
- **`granularPhaseModel`**:
  - Compressible granular density via `thermo.updateRho(p_fluid)`, limited by $\alpha \le \alpha_\text{max}$.
  - Pressure-gradient coupling term $\alpha_s\sum_f\alpha_f\nabla p_f/(1 - \alpha_s)$ (upstream: $\alpha_s\sum_f\nabla p_f$). (2026-10-12) Replaced by $\alpha_s\nabla p$ with the mixture fluid pressure: the individual fluid-phase pressures are not defined where those phases are absent, and their gradients made the Al at the unreacted charge/air interface compact and move (Section 7).
  - Blended velocity/energy/$\Theta$ decode (`0a181853` relaxed the thresholds to `residualAlpha`). (2026-10-12) Made conservative as for `fluidPhaseModel`. Where $\alpha$ reaches `alphaMax`, the grain density becomes `alphaRho/alphaMax` instead of the mass being reduced. The packing limit itself must come from the frictional stress (small `alphaDeltaMin`, Section 6.2).
- **`phaseModel` / `phasePair`**: `d0()`/`celld0()` (initial diameter) accessors; guards on $\nu$ and $d$ in the Reynolds number.

---

## 4. Thermodynamics, Chemistry & Numerics Refinements

### 4.1 Mass Fraction Handling & Mixture Robustness
- **`basicSpecieBlastMixture::correctMassFractions()`**: Cells and faces where $\sum Y_i \approx 0$ are set to the `defaultSpecie` (or species 0) instead of triggering a fatal error. `defaultSpecie` / `inertSpecie` is now also looked up when phases are named.

### 4.2 Chemistry Solvers & Combustion Models
- `blastChemistryModel` and the chemistry solvers (`blastChemistrySolver`, `noBlastChemistrySolver`, `odeBlastChemistrySolver`, `EulerImplicitBlastChemistrySolver`) take a `multicomponentBlastThermo&` and are instantiated for `phaseFluidBlastThermo` (`forCoeffGases`).
- **Selection names changed**: chemistry solver type names now include `Thermo::typeName` (`makeBlastChemistrySolver.H`).
- `reactingFluidBlastThermo` and `reactingSolidBlastThermo` hold a stored, `AUTO_WRITE` `Qdot_` field and call `correctMassFractions()` in the constructor.

### 4.3 Activation (Programmed Burn) Models
- **`nodeProgrammed`** (new, `src/thermodynamicModels/activationModels/nodeProgrammed/`): node-based programmed burn with $\lambda = (t - t_\text{start})/(t_\text{end} - t_\text{start})$ and `PROGRAMMED` / `PROGRAMMEDBETA` modes. Used by the `blastEulerFoam` `ChargeCompBAl` tutorials.
- **`programmedIgnition`**: adds `pName`, `pMin` (default $1.5\times10^5$ Pa) and `residualAlpha` (default $10^{-4}$). $\lambda = 0$ where the phase $\alpha \le$ `residualAlpha`. The `BETA` mode also requires $p >$ `pMin`; `PROGRAMMEDBETA` logic was restructured.
- **`activationModel`**: a detonation point is skipped until `time >= delay`.

### 4.4 Equation of State & Caloric Model Updates
- **`JWL`**: Maximum pressure `pMax` (default $10^{15}\,\text{Pa}$) in all constructors and operators. With `limit` on, $p$ is clamped to $[0, p_\text{max}]$.
- **`Murnaghan`**:
  - `pMax_` (default $10^{15}\,\text{Pa}$) in all constructors and operators.
  - `rhoMin` corrected to $\rho_0\max(1 - p_\text{ref}\kappa n, 10^{-4})^{1/n}$.
  - `rho(p, T)` implemented; `E()` uses the limited density.
  - (The constructor already accepted either `kappa` or `K0`; `2971df46` switched the case files from the misused `kappa` to `K0`.)
- **`AbelNobel` / `AbelNobelSG`**: covolume mixing rule changed from mole-fraction to mass-fraction weighting.
- **`perfectGas`**: Refactor only. $p$ is still bounded below by 0 when limiting, with **no upper bound**.
- **`janafBlastThermo`**: $C_p$, $H_a$, $C_{pt}$, $H_{at}$, $S$, $G_\text{std}$ and $dC_p/dT$ are evaluated at temperature-clamped values (`limit(T, T)`).
- **`blastThermo`**: New `THigh` upper temperature limit (default $10^5$ K, read when `limitT`; set it to e.g. 5000 K for JANAF gases) and `CpRef`/`CvRef`/`kappaRef`, with matching virtual interfaces in `phaseBlastThermo`, `mixtureBlastThermo` and `multicomponentBlastThermo`.
- **Temperature solve**:
  - `basicFluidBlastThermo`, `basicSolidBlastThermo` and `detonatingSolidBlastThermo` check the energy floor at `TLow` before the Newton solve, and cap $T$ at 10000 K.
  - `thermoModelI.H` detects a 2-cycle oscillation in the Newton iteration and averages it.
- **`detonatingFluidBlastThermo`**:
  - Mixed-zone temperature is a linearised $T_\text{low} + (e - e_\text{test})/C_{v,\text{mix}}$, clamped to $[T_\text{low}, T_\text{high}]$.
  - Internal-cell sound speed is $\min(c_1, c_2)$; boundary sound speed is capped at 20000 m/s.
  - Bug fix: `t1.Gamma*x2` → `t2.Gamma*x2`.
- **`twoPhaseFluidBlastThermo`**:
  - `XiSum` clamped to $\ge 10^{-6}$ in the mixture pressure.
  - Sound speed uses $\max(c^2\rho\Xi_\text{sum}, \text{small})$.
  - Per-cell `Pout` diagnostic for $|p| > 10^{11}$ or NaN, active only with the `twoPhaseFluid` debug switch (since 2026-10-10).

### 4.5 Diameter Models
- `diameterModel` and `reactingDiameterModel`:
  - New `d0()`/`celld0()` (initial diameter).
  - The reacting diameter is bounded to $[10^{-10}, d_0]$ and reset to $d_0$ where $\alpha \le 10^{-4}$.
  - (2026-10-13) $d_0$ is a registered field `d0.<phase>` (read if present, written). It was a copy of `d.<phase>` under the same name, so it was neither mapped on mesh refinement nor restored on restart; `tutorials/blastEulerFoam/reactingParticles` (AMR) crashed with a division by zero in `RanzMarshall::Ki` at 2.6e-7 s.
  - Number density uses $\pi d_0^3/6$.

### 4.6 Flux Schemes
- **`HLLC` phase flux scheme** (`src/finiteVolume/fluxSchemes/phaseFluxSchemes/HLLC/`):
  - Rewritten to the standard HLLC star-state flux $F + S(U^* - U)$ with separate owner/neighbour star pressures.
  - Explicit $S_\text{nei} \le 0$ branch, guarded $S^*$, and central-average fallback.
- `HLL`, `HLLC` and `Kurganov` phase flux schemes are also registered in the `solid` selection table.
- **`fluxScheme`**: Overloaded `NewSingle(phi, type)` and `NewMulti(phi, type)` factory methods for selecting a scheme by name. Not yet called anywhere.

### 4.7 Function Objects (`124926e3`)
- **`temporalMax`** (aliases `temporalFieldMax`, `coarseMeshMax`, `coarseFieldMax`, `coarseMax`) and **`temporalIntegral`** (aliases `timeIntegral`, `coarseMeshIntegral`) in `src/functionObjects/`:
  - Per-cell temporal maximum and time integral (e.g. impulse, with `pRef` and `positiveOnly`) over an optional time window.
  - Results are mapped back to coarse-mesh ancestor cells under AMR, and radial profiles are written.

---

## 5. Build System
- **`src/finiteVolume`**:
  - `Make/options` adds include paths for `physicalProperties`, `thermophysicalModels/{basic,specie,multicomponentThermo}`, blast `numerics` and blast `thermodynamicModels`.
  - It links `physicalProperties`, `fluidThermophysicalModels`, `specie` and `blastNumerics`.
  - These were added for `pressurePhaseChange`. Since `c3294dc1` no `src/finiteVolume` source uses them; the stale `reactionThermo` include (absent in OpenFOAM-12) was removed.
- **`src/thermodynamicModels/Make/files`**: compiles `multicomponentPhaseFluidBlastThermos.C`, `reactingPhaseFluidBlastThermos.C`, `activationModels/nodeProgrammed/nodeProgrammedActivation.C` and (since `c3294dc1`) `fvModels/pressurePhaseChange/pressurePhaseChange.C`.
- **Dependency rule**: `libblastThermodynamics` links against `libblastFiniteVolume` and is built after it (`src/Allwmake`). Code in `src/finiteVolume` must therefore not use thermodynamicModels classes. `pressurePhaseChange` did, which left undefined `multiphaseFluidBlastThermo` symbols and broke the link of utilities such as `addEmptyPatch` (fixed in `c3294dc1`).

---

## 6. Tutorials, Verification Cases & Tooling

### 6.1 Tutorials
New or changed since `71b97fd`:
- **`tutorials/blastFoam/pressurePhaseChangeTest`**: Minimal `pressurePhaseChange` demonstration (air to N2 bubble, `pThreshold 2e5`, `pRange 2e4`).
- **`tutorials/blastFoam/ChargeCompBAl_1D` / `ChargeCompBAl_1D_AlInert`**: 1D aluminised Composition B detonation and afterburning.
- **`tutorials/blastEulerFoam/ChargeCompBAl`** (2D axisymmetric wedge), **`ChargeCompBAl_1D`**, **`StationaryCharge_1D`**, **`axisymmetricCharge`**: multiphase Euler–Euler aluminised charges (`nodeProgrammed` activation).
- **`tutorials/blastFoam/multicomponent_shock_bubble`**: Shock interaction with a multicomponent gas bubble.
- **`tutorials/blastFoam/multiphase_reacting_bubble` / `multiphase_reacting_bubble_1d`**: Multiphase reacting bubble ignition and expansion.
- **`tutorials/blastFoam/axisymmetricChargeMod` / `axisymmetricChargeReacting`**: Axisymmetric detonation with reacting multicomponent products.
- **`tutorials/blastFoam/shock_Air_Al_Air`**: Two-phase (`air`/`aluminum`) air–Al–air shock tube.
- **`tutorials/blastReactingFoam/obliqueShock`**: Added `Allrun.parallel`, `decomposeParDict` and `dynamicMeshDict`.
- Minor: `tutorials/blastFoam/building3DWorkshop/constant`, `tutorials/blastReactingFoam/reacting/constant/dynamicMeshDict`.

### 6.2 Verification Cases (`tests/`)
- **0D**:
  - `phaseTransfer_0D` (isochoric eJWLC products → JANAF gas) and `phaseTransfer_0D_3phase` (`2971df46`: adds inert Al). Since 2026-10-10 both run **blastEulerFoam** (replacing the blastFoam versions): products as a `basic` eJWLC phase, gas as a reacting JANAF/Abel–Nobel phase, Al as a granular Murnaghan phase with RanzMarshall heat transfer. Both are in Tier 1 of `bin/run_blast_tests.py`.
  - `setup_case.py [--three-phase]` regenerates the per-phase `0/` fields (all phases at 100 bar and `T0 = 825.76` K); `verify_0d_phase_transfer.py [--three-phase]` checks mass, momentum and total energy (sum over phases, tolerance 1e-8), pressure (5%) and conversion.
  - `verify_al_particle_kinetics.py` and `verify_highT_air_thermo.py` are standalone numpy models (they do not run a solver).
- **1D**:
  - `multiphase_advection`, `planar_znd`, `shock_dissociating_air` (`4cb49946`); `planar_znd_phaseTransfer` (`0a181853`); `gunib1d_houim1/2/3` (`124926e3`).
  - Each has a `verify_*.py`. The planar ZND verifiers need the EuDetCode build in `external/eudetcode` (a separate checkout), and `verify_gunib1d_benchmarks.py` a local gunIB1D build.
- **1D_spherical**:
  - `spherical_blast` (C4, Kingery-Bulmash, `4cb49946`).
  - `ChargeCompBAl_1D_spherical` (`0a181853`) and `_reacting` (`124926e3`).
  - `_alCombustion` (`6f693f7b`; gas-phase Al reactions Al+O2/CO2/H2O and H2+O2 added in `4fd73e06`).
  - (2026-10-12) `_interface` (`verify_spherical_interface.py`, Tier 3): the spherical charge never initiated, for 5 µs. Its exact solution is the initial state. It checks that the unreacted charge/air interface stays at rest and conserves mass and energy.
- **Case settings changed in `2971df46`**:
  - Murnaghan `kappa` → `K0`.
  - compBAl `residualRho 1e-6`.
  - Gas `THigh 5000`.
  - eJWLC products `Cv 153.3` → `1450`.
  - Charge densities made consistent with reactant $\rho_0$.
- `68cbb94b` and `bdccfc0e` only extend run times.
- **`alphaDeltaMin`** (2026-10-13): with the upstream tutorial value 0.05 the JohnsonJackson frictional pressure is capped at about 10 kPa and cannot stop compaction. `gunib1d_houim3` now uses $10^{-6}$ (compaction speed 183.8 m/s vs 185.25 m/s reference; the bed previously packed to `alphaMax`, and the original decode removed 41% of the particle mass). The Comp-B/Al tests (`1D_spherical/ChargeCompBAl_1D_spherical*`) and the 1D tutorials (`ChargeCompBAl_1D`, `StationaryCharge_1D`) use $10^{-3}$; the Al previously packed to `alphaMax` at the charge centre. The 2D `ChargeCompBAl` and the other upstream tutorials are unchanged.
- **Products/gas pressure and temperature relaxation** (2026-10-13): the inert and reacting spherical cases, `_interface` and the 1D tutorials set `alphaDilute 1` and a `(compBAl in gas)` `constantCoefficient` heat transfer with `Kh 1e15` (instantaneous, through the exact exponential update in `phaseSystem::relaxTemperature()`), so numerically mixed contact cells are in pressure and temperature equilibrium. `_alCombustion` keeps the dilute-only setting (its Al burn puts trace gas inside the charge near its covolume limit; the stronger relaxation crashes it at breakout). `verify_spherical_interface.py` checks the air pressure on cells without explosive or Al.

### 6.3 Tooling
- **`bin/run_blast_tests.py`**: tiered regression runner (`--tier 0|1|2|3|all`, `--case`, `--report`). It writes `blast_test_report.md`. Paths are hard-coded to `/opt/OpenFOAM/...`.
- **`utilities/thermo_fitter/fit_ejwlc_nobleabel.py`**: eJWLC / Abel-Nobel thermodynamic reconciliation fitter.
- **`utilities/thermo_fitter/compute_omegaT.py`** (2026-10-10): computes and optionally writes (`--write`) the eJWLC `omegaT` for a case: solves the products isentrope for $V_\text{trans}$ at the transfer pressure, takes $T_\text{trans}$ from the gas EOS with the transfer composition (`--gas-case` to take the gas from another case, `--inert Y:RHO` for fits that include inert material), and $\omega_T = \ln(T_\text{CJ}/T_\text{trans})/\ln(V_\text{trans}/V_\text{CJ})$.
- **`utilities/generate_report_figures.py`**: figures for `report/TBX-blastFoam.tex`.

---

## 7. Known Issues & Fixes

### Fixed (2026-10-10)
- **`isMoleFraction` was ignored** (`pressurePhaseChange.C` discarded the lookup result, so mole fractions were silently applied as mass fractions). Mole fractions are now read and converted to mass fractions (Section 2.1). Verified: the 0D case run with equivalent mole fractions reproduces the mass-fraction run exactly. Affected inputs: `tutorials/blastEulerFoam/{ChargeCompBAl,ChargeCompBAl_1D,StationaryCharge_1D}/constant/fvModels.disabled`.
- **`reactingDiameterModel` `d0`** (2026-10-13): not mapped on mesh refinement, which crashed `tutorials/blastEulerFoam/reactingParticles` (Section 4.5).
- **Gas `residualRho` floor** (2026-10-13): removed outside residual cells (Section 3.3).
- **`pressurePhaseChange` `lambda` factor** (2026-10-13): the `addSup(rho, field, eqn)` energy/alpha path (used only without a time integrator) omitted the source-phase `lambda`. All paths now use one `mDot()` helper; the time-integrator results are unchanged.
- **Debug output**: the per-cell `CELL DEBUG` print in `multiphaseFluidBlastThermo::calculate()` was removed. The extreme-pressure print in `twoPhaseFluidBlastThermo` is now behind its debug switch.
- **Energy offset**: `energyOffset()` was meant to compute a state-dependent offset $\sum_k f_k h_{e,k}(T) - h_{e,\text{source}}(T)$ at the shared temperature. It never ran, because it looked the thermo up under the wrong registry name (`"thermophysicalProperties"` instead of `physicalProperties`).
  - When enabled for testing, it added energy in proportion to the JWL temperature error: about 480 kJ/kg and +56% pressure in the 0D case with the power-law isentrope temperature. A state-dependent offset is not a datum shift and is not energy-conservative.
  - It was replaced by a constant datum shift only (`eOffset`, default 0; Section 2.1). Results of existing cases are unchanged, since the dynamic branch never ran.
  - The datums are reconciled through the products `E0` (report Section 3.3).
- **`tests/0D/phaseTransfer_0D` pressure check** (16.48% after `2971df46`'s products `Cv 153.3 → 1450`):
  - Cause (blastFoam version): the single mixture energy shares one temperature between the phases, and the eJWLC power-law isentrope put the products at 1318.69 K where the gas EOS gives about 826 K. As the shared temperature relaxed, energy moved between the phases in proportion to $C_v$.
  - Temperature fix: new optional eJWLC exponent `omegaT` (Section 1.2), keeping the EuDetCode `T_CJ`. `omegaT` was computed with `compute_omegaT.py` and set in all eJWLC cases: 0.355493 (products-only Comp-B/Al fit: 0D tests, `spherical_alCombustion`, `blastFoam/ChargeCompBAl_1D`), 0.375928/0.375955 (full Comp-B/Al fit with 27% inert Al: `spherical`, `spherical_reacting`, blastEulerFoam tutorials), 0.469442 (Comp-B into air: `planar_znd`, `planar_znd_phaseTransfer`), 0.353900 (`ChargeCompBAl_1D_AlInert`).
  - The 0D tests now run blastEulerFoam (independent phase energies, so the shared-temperature effect no longer applies): 0.58% maximum pressure deviation, gas formed at 830.6 K, momentum exactly zero, energy conserved to about 5e-9 (since the conservative decode; before it, the dilute-phase reset at 91 µs removed 1.8% / 1.1%).
- **`pressurePhaseChange` in blastEulerFoam moved mass without momentum or energy**: the source kept its momentum and energy while losing mass, and the targets received mass with neither. Now the transferred mass carries the source velocity and total energy (Section 2.1).
- **`tests/0D/phaseTransfer_0D_3phase`** (107.9 bar initial pressure in the blastFoam version, from the mixture-energy closure with an energy-independent Murnaghan phase): in blastEulerFoam the phases have their own energies; the case starts at 100 bar and its pressure check passes (0.58%).
- **Dangling `tmp` references** in `RanzMarshallHeatTransfer.C` and `sphericalHeatTransfer.C` (bound a reference to the temporary returned by `d()`).
- **EuDetCode paths**: the test scripts and the fitter pointed at a non-existent `modules/eudetcode`. They now use `external/eudetcode`, a separate git checkout not tracked by this repository. `tests/1D/verify_planar_znd.py` passes against it (detonation velocity within 0.076%). Tier 0 still needs the EuDetCode GoogleTest binary (`build/tests/eudetcode_tests`), which is not built.
- **0D verifier** (`verify_0d_phase_transfer.py`): energy conservation is checked on the datum-corrected energy, stored energy minus `eOffset` × transferred mass.

- **Non-conservative `fluidPhaseModel` / `granularPhaseModel` decode** (2026-10-12): dilute phases ($\alpha < 0.01$) had their density, velocity and energy reset with the conserved variables rewritten. This destroyed 1.8% / 1.1% of the energy in the 0D tests and deleted the energy of the products expanding into the air in the spherical cases. Removing the resets exposed two causes of spurious GPa pressures that the resets had masked:
  - the granular buoyancy term built from the individual fluid-phase pressures (Section 3.3);
  - trace fluid phases at fluid/fluid contacts whose decoded density drifts (isostrain without pressure relaxation; fixed by `equilibrateDiluteVolumes()`).
  Both decodes are now conservative. 0D energy is conserved to ~5e-9, and the spherical verifiers (0.64% each since 2026-10-13), the new interface test and the blastEulerFoam tutorials pass. With the gas density floor also removed (2026-10-13), the spherical air shock reaches 10 cm at 17.9 µs with 14 MPa (originally 25.4 µs and 90 MPa). Details, rejected approaches and remaining points: `report/dilute_phase_conservation_plan.md`.

### Open
- **Pressure non-equilibrium in non-dilute mixed fluid cells**: without pressure relaxation, fluid phases with $\alpha > 0.01$ in the same cell undergo the same volumetric strain. At the products/air contact about 50 cells are mixed, with the products at 2–5000 times the gas pressure and the entrained air rarefied at its 5000 K ceiling (2026-10-13). The upstream `instant` pressure relaxation collapses the time step; extending `equilibrateDiluteVolumes()` to all phases (`alphaDilute 1`) equalises the pressures but leaves trace-mass products at 1e5 K. Adding products/gas temperature relaxation fixes this in the inert and reacting cases (Section 6.2); still open for `_alCombustion`, and the common temperature in the planar tutorial's mixed cells (about 9000 K) is high. See `report/dilute_phase_conservation_plan.md` (Section 6).
- **`tutorials/blastEulerFoam/reactingParticles`**: after the `d0` fix it runs past 2.6e-7 s, but the time step becomes erratic within the first few hundred steps (down to 1e-15 s) and collapses to 1e-22 s near 2.7e-4 s (of 2.5e-3 s), followed by a floating-point exception in the gas thermo. The original decode behaves the same (collapse at 2.65e-4 s). Not analysed.
- **Al vapour inside the charge (`_alCombustion`)**: the burning Al releases its mass as Al vapour into the gas (air) phase, also inside the detonating charge, where that phase is otherwise absent. The vapour ends up at 4500–8000 kg/m³ and 5–20 GPa, outside the Abel–Nobel/JANAF range; this is why the case cannot use the products/gas pressure–temperature relaxation. Needs a modelling decision (see `report/dilute_phase_conservation_plan.md`, Section 3).
- **gunIB1D comparison**: the local gunIB1D build (`/home/kevnor/Codes/gunIB1D/gunIB1D`) stops with a floating-point exception on all three Houim cases, so the WP 2.3 cross-solver comparison cannot be repeated. (2026-10-13) `verify_gunib1d_benchmarks.py` now finds the binary (`GUNIB1D_DIR`, `build/gunIB1D` or `gunIB1D`) and the `Validation/Houim*.info` cases, runs gunIB1D in a temporary directory instead of the gunIB1D checkout, and first checks Houim 3 against the Houim (2016) reference without gunIB1D (compaction speed within 2%, particle mass to 1e-8; currently 0.81% and 4e-12). It still fails while gunIB1D crashes.
- **`tutorials/blastFoam/ChargeCompBAl_1D_AlInert`** sets `E0 = -5.484e6` on the products and `eOffset = 5.484e6` in `fvModels`. The two cancel, so the transferred energy is as if `E0 = 0`: the datums are not reconciled for its `hConst` gas. It should use a calibrated `E0` and no `eOffset`.
- **Possible double counting of Al** in the full Comp-B/Al cases (`spherical`, `spherical_reacting`, blastEulerFoam `ChargeCompBAl*`, `StationaryCharge_1D`): their eJWLC fit (ρ0 = 1820) includes the 27% inert Al in its mass, while they also carry a separate Al particle phase.
- **blastFoam mixture closure**: in blastFoam two-/multiphase systems each phase EOS is evaluated with the mixture energy (Allaire-type closure, exact for Mie–Grüneisen phases only); a phase whose pressure is independent of energy (Murnaghan) is then inconsistent. Not relevant to the blastEulerFoam problems.
