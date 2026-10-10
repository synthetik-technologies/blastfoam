# Library Updates (Compared to `origin/develop` @ 71b97fd)

This document summarises the updates, new features, thermodynamic models, numerical changes and bug fixes in this fork relative to `origin/develop` (commit `71b97fd7c9195fdb342ed4230d258dd1dde80b55`, the last merge of upstream `synthetik-technologies/blastfoam`).

*Last synchronised with the code at commit `c3294dc1` (2026-10-10).* Commits covered: `4cb49946`, `0a181853`, `124926e3`, `6f693f7b`, `4fd73e06`, `2971df46`, `68cbb94b`, `bdccfc0e`, `8b274c09`, `fa2c2283`, `c3294dc1`. See [Known issues](#7-known-issues) for open defects found while reviewing these changes.

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
- `targetFractions` sets the target-phase species distribution. **These are always applied as mass fractions**; see Known issues for `isMoleFraction`.
- Inter-phase energy offset $e_\text{off} = \sum_i Y_i\,h_{e,i}^\text{target}(T) - h_e^\text{source}(T)$, evaluated per cell and per boundary face, plus an optional constant `eOffset`.
- Adds source terms to the phase continuity equations ($\alpha_i$, $\rho_i$, $\alpha\rho_i$) and to the energy equation:
  - blastFoam (single mixture energy): net `rhoE` source $-\dot m\,e_\text{off}$.
  - blastEulerFoam (`124926e3`): per-phase `alphaRhoE.<phase>` sources $\pm\dot m\,e_\text{off}$.
  - No momentum source (`addSup` for vectors is a no-op).
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
- **Heat transfer** (`RanzMarshall`, `spherical`): $K$ scales as $d/\max(d, d_0)^3$ instead of $1/d^2$.
- **Kinetic theory**:
  - `kineticTheorySystem` collisional `lambda` term rewritten per cell with residual-$\alpha$ guards.
  - Bounded denominators in the `CarnahanStarling` and `SinclairJackson` radial models.
  - $\sqrt{\max(\Theta, \text{small})}$ in Gidaspow/Syamlal conductivity and viscosity.
  - **Syamlal viscosity: the collisional term expression was changed** (a formula change, not only a guard).
- **`phaseSystem`**:
  - Velocity, pressure–temperature-equilibrium and temperature relaxation exchange are zeroed where either phase is residual.
  - Mixture temperature averages only active phases (fallback 288.15 K); mixture density is limited.
  - `decode()` decodes granular phases first, rescales the carrier volume fraction, then decodes fluid phases.
- **`fluidPhaseModel` decode** (note: these modify the state non-conservatively):
  - Density is blended towards `residualRho` for $\alpha < 10^{-3}$ (adds mass).
  - Velocity is capped at 20000 m/s, and zeroed at small $\alpha$.
  - Energy is clamped between $h_e(\max(100\,\text{K}, T_\text{low}))$ and $h_e(T_\text{high})$.
  - Energy and temperature are reset to 288.15 K where $\alpha \le 0.01$.
  - `solveImplicit` calls `decode()` instead of `thermo.correct()`.
- **`granularPhaseModel`**:
  - Compressible granular density via `thermo.updateRho(p_fluid)`, limited by $\alpha \le \alpha_\text{max}$.
  - Pressure-gradient coupling term $\alpha_s\alpha_f\nabla p/(1 - \alpha_s)$.
  - Same blended velocity/energy/$\Theta$ decode. (`0a181853` relaxed the thresholds to `residualAlpha`.)
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
  - Prints a `Pout` diagnostic for every cell with $|p| > 10^{11}$ or NaN (not throttled).

### 4.5 Diameter Models
- `diameterModel` and `reactingDiameterModel`:
  - New `d0()`/`celld0()` (initial diameter).
  - The reacting diameter is bounded to $[10^{-10}, d_0]$ and reset to $d_0$ where $\alpha \le 10^{-4}$.
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
  - `phaseTransfer_0D` (isochoric eJWLC products → JANAF gas).
  - `phaseTransfer_0D_3phase` (`2971df46`: adds inert Murnaghan Al).
  - `setup_case.py [--three-phase]` regenerates the `0/` fields; `verify_0d_phase_transfer.py [--three-phase]` checks the results.
  - `verify_al_particle_kinetics.py` and `verify_highT_air_thermo.py` are standalone numpy models (they do not run a solver).
- **1D**:
  - `multiphase_advection`, `planar_znd`, `shock_dissociating_air` (`4cb49946`); `planar_znd_phaseTransfer` (`0a181853`); `gunib1d_houim1/2/3` (`124926e3`).
  - Each has a `verify_*.py`. The planar ZND verifiers and `verify_gunib1d_benchmarks.py` need the external EuDetCode and gunIB1D builds.
- **1D_spherical**:
  - `spherical_blast` (C4, Kingery-Bulmash, `4cb49946`).
  - `ChargeCompBAl_1D_spherical` (`0a181853`) and `_reacting` (`124926e3`).
  - `_alCombustion` (`6f693f7b`; gas-phase Al reactions Al+O2/CO2/H2O and H2+O2 added in `4fd73e06`).
- **Case settings changed in `2971df46`**:
  - Murnaghan `kappa` → `K0`.
  - compBAl `residualRho 1e-6`.
  - Gas `THigh 5000`.
  - eJWLC products `Cv 153.3` → `1450`.
  - Charge densities made consistent with reactant $\rho_0$.
- `68cbb94b` and `bdccfc0e` only extend run times.

### 6.3 Tooling
- **`bin/run_blast_tests.py`**: tiered regression runner (`--tier 0|1|2|3|all`, `--case`, `--report`). It writes `blast_test_report.md`. Paths are hard-coded to `/opt/OpenFOAM/...`.
- **`utilities/thermo_fitter/fit_ejwlc_nobleabel.py`**: eJWLC / Abel-Nobel thermodynamic reconciliation fitter.
- **`utilities/generate_report_figures.py`**: figures for `report/TBX-blastFoam.tex`.

---

## 7. Known Issues
Found while synchronising this document (2026-10-10); not yet fixed.
- **`pressurePhaseChange` ignores `isMoleFraction`**: `pressurePhaseChange.C:662` discards the result of `lookupOrDefault`, so `isMoleFraction_` stays `false` and `targetFractions` are always treated as mass fractions. Cases that set `isMoleFraction true` silently get the wrong composition (`tutorials/blastEulerFoam/{ChargeCompBAl,ChargeCompBAl_1D,StationaryCharge_1D}/constant/fvModels.disabled`).
- **`pressurePhaseChange` does not conserve total energy exactly**:
  - The blastFoam path adds a net energy source $-\dot m\,e_\text{off}$.
  - The per-phase path moves only $\pm\dot m\,e_\text{off}$, not the transferred mass's own internal and kinetic energy.
  - The `lambda.<phase>` factor is applied to the mass sources but not in the `addSup(rho, field, eqn)` (energy/alpha) path.
- **Debug output left in**: `multiphaseFluidBlastThermo.C:63` prints a `CELL DEBUG` line for every cell with $T < 200$ K, and `twoPhaseFluidBlastThermo` prints every extreme-pressure cell.
- **`tests/0D/phaseTransfer_0D` fails its pressure check** since `2971df46`: 16.48% deviation vs the 5% limit (passes at 4.09% with the old `Cv 153.3`).
