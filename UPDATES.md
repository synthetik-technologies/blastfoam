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
  - (2026-10-13) Instantaneous velocity relaxation (`VRelaxation instant`; both upstream call sites now use `relaxFluidVelocitiesInstant()`). The per-stage block in `solve()` added $\tfrac12\alpha\rho|V_I - U|^2$ to the total energy, which changes the internal energy by $U\cdot(U - V_I)$: it created the dissipated energy a second time and heated (or, for $0 < U < V_I$, cooled) each phase by its own velocity defect. The block in `relaxVelocity()` kept the internal energies and so destroyed the dissipation. Now all fluid phases take the mass-averaged velocity and each keeps its internal energy plus its own dissipation $\tfrac12|U - V_I|^2$, which conserves total energy. A phase with little mass can carry momentum far beyond its energy (its face fluxes are weighted by the face volume fraction, Section 7); where its internal energy would fall below the thermo range, that kinetic energy is not counted, and heating above the upper thermo limit is not applied. Sharing the dissipation by mass instead was tried: it keeps the stationary interface stable for 40 µs, but moves the spurious kinetic energy of the gas in mixed product cells into the products (70 000 K, the reacting case fails at 12 µs).
  - (2026-10-13) Drag work: the work of the implicit drag momentum exchange is evaluated at the interfacial velocity averaged over the update (old and post-drag velocities). With the interfacial velocity of the particles (`single`), the gas lost $\Delta M\cdot U_{p,\text{old}}$ while the particle kinetic energy rose by $\Delta M\cdot\bar U_p$, creating $|\Delta M|^2/(2\alpha\rho_p)$ per step.
- **`fluidPhaseModel` decode**:
  - (`4cb49946`) Density blending towards `residualRho` below $\alpha = 10^{-3}$, velocity blending to zero below $10^{-2}$, and an energy/temperature reset to 288.15 K below $\alpha = 0.01$, all written back to the conserved variables (non-conservative). `solveImplicit` calls `decode()` instead of `thermo.correct()`.
  - (2026-10-12) Made conservative: `alphaRho`, `alphaRhoU` and `alphaRhoE` are only modified where the phase is residual ($\alpha \le$ `residualAlpha`, reference state 288.15 K), as upstream. The velocity cap (20000 m/s) and the energy clamp to the thermo temperature range act on the primitive variables only.
  - (2026-10-13) The decoded density is no longer bounded below by `residualRho` (upstream `rho.max(residualRho)`); `residualRho` is only the reference density of residual cells. The gas `residualRho 1.2` held expanded gas at ambient density in up to 90 cells of the spherical runs and overstated its pressure up to twentyfold.
- **`granularPhaseModel`**:
  - Compressible granular density via `thermo.updateRho(p_fluid)`, limited by $\alpha \le \alpha_\text{max}$.
  - Pressure-gradient coupling term $\alpha_s\sum_f\alpha_f\nabla p_f/(1 - \alpha_s)$ (upstream: $\alpha_s\sum_f\nabla p_f$). (2026-10-12) Replaced by $\alpha_s\nabla p$ with the mixture fluid pressure: the individual fluid-phase pressures are not defined where those phases are absent, and their gradients made the Al at the unreacted charge/air interface compact and move (Section 7).
  - Blended velocity/energy/$\Theta$ decode (`0a181853` relaxed the thresholds to `residualAlpha`). (2026-10-12) Made conservative as for `fluidPhaseModel`. Where $\alpha$ reaches `alphaMax`, the grain density becomes `alphaRho/alphaMax` instead of the mass being reduced. The packing limit itself must come from the frictional stress (small `alphaDeltaMin`, Section 6.2).
  - (2026-10-13) Thermal energy transport. `alphaRhoE` is the grain internal energy, but the HLLC flux scheme supplied the total-energy flux $(\rho(e + \tfrac12|U|^2) + p)U$ (upstream as well), so the kinetic energy was advected into the thermal energy, and the grain compression work was missing although the grain density follows the mixture pressure (Murnaghan cold energy $E_c(\rho)$, 351 kJ/kg at 15 GPa). The Al cooled behind the detonation to the 100 K clip. Now `alphaRhoE` holds the thermal part $\alpha\rho(e - E_c(\rho))$, advected with the mass flux (reconstructed `e` upwinded by `alphaRhoPhi`), and decode adds $E_c$ at the new density (`eCold()`, evaluated through the thermo at the reference density from `updateRho()`; `blastThermo::rho0()` is the old-time density). The Al stays at 288 K under compression and is then heated by the products (316 K at 10 µs). The gunIB1D comparison is unchanged (constant-density particles).
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
  - (2026-10-13) Optional oxidiser limit (`oxidiserPhase`, `oxidiserSpecies`, `YoxidiserRef`): the burn rate is multiplied by the fraction of the surrounding fluid volume occupied by the oxidiser phase, $\alpha_\text{ox}/(1-\alpha_d)$, and by $\min(1, \sum Y_\text{ox}/Y_\text{ref})$. Without `oxidiserPhase` the rate is unchanged. Used by `_alCombustion` (gas phase; O2, CO2, H2O; $Y_\text{ref} = 0.24$).

### 4.6 Flux Schemes
- **`HLLC` phase flux scheme** (`src/finiteVolume/fluxSchemes/phaseFluxSchemes/HLLC/`):
  - Rewritten to the standard HLLC star-state flux $F + S(U^* - U)$ with separate owner/neighbour star pressures.
  - Explicit $S_\text{nei} \le 0$ branch, guarded $S^*$, and central-average fallback.
  - A phase absent on one side of a face ($\alpha \le 10^{-6}$) is treated as vacuum (upstream). (2026-10-13) Tried: giving the absent side the state of the present side. This removes the start-up transient of the stationary charge/air interface (pressure stays at 1.01325 bar, velocities below 0.1 m/s for 40 µs), but at a detonating contact the products' face pressure flux then carries their full pressure work into the air cell while the balancing $P_I\nabla\alpha$ term there uses the air pressure, and the trace products reach 1e5 K (the spherical cases fail at 9–11 µs). Not adopted; see Section 7.
  - (2026-10-14) **Trace sides take the mixture impedance** (`alphaTrace`, `phaseProperties`, default 0: off). Where a fluid phase is a trace on one side of a face ($\alpha < $ `alphaTrace`, blended linearly in $\alpha$/`alphaTrace`; at $\alpha \le 10^{-6}$ in place of the vacuum state), the wave speeds and a common star pressure of that phase's Riemann problem use the cell's mixture density, frozen sound speed $\sqrt{\sum\alpha_k\rho_k c_k^2/\sum\alpha_k\rho_k}$, velocity and pressure on that side; the star flux is written directly as $\rho^* S^* U^* + p^* n$ and $(\rho E^* + p^*)S^*$. The phase's own density and energy still define the fluxes, so they stay conservative; where the phase is above `alphaTrace` on both sides the flux is unchanged. Before, a trace of unreacted explosive in the first air cell had the explosive's impedance (4.5e6 Pa·s/m): when the air there moved outwards by 0.015 m/s, the explosive face pressure dropped by 0.35 bar, and with the upwind face volume fraction switching from the trace to 0.818 this pulled the air cell back by 0.6 m/s in one step (the interface instability, Section 7). The mixture state comes from `phaseSystem` through `phaseFluxScheme::setMixture()`; only `HLLC` uses it.
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
  - `setup_case.py [--three-phase]` regenerates the per-phase `0/` fields (all phases at 100 bar and `T0 = 836.10` K since 2026-10-13, 825.76 K before); `verify_0d_phase_transfer.py [--three-phase]` checks mass, momentum and total energy (sum over phases, tolerance 1e-8), pressure (5%) and conversion.
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
- **Products/gas pressure and temperature relaxation** (2026-10-13): the inert and reacting spherical cases, `_interface` and the 1D tutorials set `alphaDilute 1` and a `(compBAl in gas)` `constantCoefficient` heat transfer with `Kh 1e15` (instantaneous, through the exact exponential update in `phaseSystem::relaxTemperature()`), so numerically mixed contact cells are in pressure and temperature equilibrium. `_alCombustion` kept the dilute-only setting until 2026-10-14: with the stronger relaxation the unreacted charge/air interface became unstable at about 3 µs and put a 0.8 kbar wave ahead of the detonation (measured speed +9.9%). Since the interface fix (`alphaTrace`) it uses the relaxation too (see below). `verify_spherical_interface.py` checks the air pressure on cells without explosive or Al.
- **Products-only fit and Al volume fraction** (2026-10-13): every case with a separate Al phase now uses the EuDetCode products-only eJWLC fit (`--jwl-products-only`; ρ0 = 1624.45, the fit `_alCombustion` already used) and the reactant Murnaghan ρ0 = 1624.45, with α_compBAl = 0.818 and α_Al = 0.182 in the charge (27% inert Al by mass at 1820 kg/m³). Before, the inert and reacting spherical cases, `_interface` and the blastEulerFoam tutorials `ChargeCompBAl`, `ChargeCompBAl_1D` and `StationaryCharge_1D` used the full fit (ρ0 = 1820, which already contains the Al) with α_Al = 0.3, i.e. 2.5 times the Al mass in total. Inert case: peak 17.91 GPa in the charge (12.74 before), 2.85 GPa at breakout (1.40), 38 MPa at 80 mm, 27 MPa at 100 mm; D = 7664 m/s (0.68%).
- **Al burns only with gas-phase oxidiser** (`_alCombustion`, 2026-10-13): the `pressureBased` rate depends on pressure alone, so the Al burned throughout the detonated charge and released vapour into the (absent) gas phase there. With the oxidiser limit (Section 4.5) there is no regression where the Al is surrounded by products (≤ 0.008% at 20 µs), and in 1D almost none in 40 µs: the Al cloud stays inside the products ball and only its α ~ 1e-4 edge meets shocked air (1.1e-4 of the Al mass burned between 10 and 40 µs). `verify_spherical_compBAl_alCombustion.py` checks this instead of the former in-charge regression (min d within 0.85–0.99 d0). Peak 18.25 GPa, breakout 2.92 GPa, D = 7661 m/s (0.64%).
- **Fits from the current EuDetCode** (2026-10-13): the products-only fit was regenerated with the updated EuDetCode (`CJDetonation "RDX:0.42 TNT:0.28 Al:0.3" 1820 BKW --inert "Al(s):0.9" --jwl-products-only`, gas database `Gas24CHNO`): D = 7154.23 m/s, p_CJ = 20.15 GPa (was 7612.3 m/s, 21.53 GPa), T_CJ = 3844.39 K, ρ_CJ = 2150.55; A = 224.943 GPa, B = 2.77303 GPa, C = 1.47024 GPa, R1 = 3.456, R2 = 0.813403, ω = 0.302717; with `products_only_fit.py` (Section 6.3) E0 = −503470.78 J/kg, omegaT = 0.371065, initiation e0 = 5.0174e6 J/kg, molWeight 24.02 (EuDetCode active molar mass). Cv stays 1450 J/(kg K) (JANAF Cv of the freeze-out gas at 2500 K; EuDetCode's BKW value is 1792). Applied to the 0D tests, the four spherical cases and the blastEulerFoam tutorials `ChargeCompBAl`, `ChargeCompBAl_1D`, `StationaryCharge_1D` (vDet 7154.23, Pcj 2.01512e10); the verifiers use the new D and p_CJ. The transfer composition (0D `fvModels`, the `fvModels.disabled` of `_reacting` and `_alCombustion`) is the new freeze-out gas, which includes H2 (0.43%; added to the 0D gas thermo); the 0D initial state is 35.233992 kg/m³ at 836.10 K. Results: 0D pressure constant to 3e-7 (0.58% before; the gas now forms at the products temperature), energy to 4.2e-9 / 2.2e-9 with Al. Spherical inert case: D = 7162 m/s (0.11%), peak 16.90 GPa, breakout 2.96 GPa, 44 MPa at 80 mm, 31 MPa at 100 mm (arrival 14.3 µs), Al to 2180 m/s; reacting: breakout 3.28 GPa; `_alCombustion`: peak 19.26 GPa, breakout 2.82 GPa, regression in products-only cells ≤ 0.012%, 1.8e-4 of the Al mass burned between 10 and 40 µs. The full-mixture fit (no `--jwl-products-only`) is A = 1428.32 GPa, B = 12.1675 GPa, C = 1.13315 GPa, R1 = 5.85811, R2 = 1.38613, ω = 0.301268 (report Table 1); no case uses it.

- **`alphaTrace 1e-3`** (2026-10-14) in all Comp-B/Al cases with a separate Al phase (the four spherical tests and the blastEulerFoam tutorials `ChargeCompBAl`, `ChargeCompBAl_1D`, `StationaryCharge_1D`; Section 4.6). `_interface`: largest phase velocity 0.0009 m/s (1.8 m/s before), pressure 101325 Pa, mass and energy changes about 1e-9; run to 40 µs at both write intervals it stays below 1.4 mm/s and 1.01 bar (before: a 114 bar event at 14 µs in one of two runs). Inert and reacting cases: peak 16.61 GPa in the charge, 3.32 GPa at breakout, 44 / 29 / 21.0 MPa at 80 / 100 / 150 mm, D = 7162.88 m/s (0.12%). `_alCombustion` (2026-10-14, now also with `alphaDilute 1` and the `Kh 1e15` products/gas heat transfer): 16.61 GPa in the charge (19.26 GPa before, raised by the interface wave), 3.33 GPa at breakout, 44 / 29 / 21.0 MPa, the same as the inert case; Al burnt by 40 µs 0.023%. With the dilute-only setting it gave 3.47 GPa and 52 / 32 / 20.3 MPa with the air shock 1 µs earlier, an effect of the products/gas pressure non-equilibrium (ratios up to 4.7e4, products at the 1e5 K cap in mixed cells), not of the Al burn (0.012%). Inert case also run with `alphaTrace` 0, 1e-4 and with maxCo 0.19: all reach 40 µs. 0D and gunIB1D results unchanged.

### 6.3 Tooling
- **`bin/run_blast_tests.py`**: tiered regression runner (`--tier 0|1|2|3|all`, `--case`, `--report`). It writes `blast_test_report.md`. Paths are hard-coded to `/opt/OpenFOAM/...`.
- **`utilities/thermo_fitter/fit_ejwlc_nobleabel.py`**: eJWLC / Abel-Nobel thermodynamic reconciliation fitter.
- **`utilities/thermo_fitter/compute_omegaT.py`** (2026-10-10): computes and optionally writes (`--write`) the eJWLC `omegaT` for a case: solves the products isentrope for $V_\text{trans}$ at the transfer pressure, takes $T_\text{trans}$ from the gas EOS with the transfer composition (`--gas-case` to take the gas from another case, `--inert Y:RHO` for fits that include inert material), and $\omega_T = \ln(T_\text{CJ}/T_\text{trans})/\ln(V_\text{trans}/V_\text{CJ})$.
- **`utilities/thermo_fitter/products_only_fit.py`** (2026-10-13): blastEulerFoam entries from an EuDetCode products-only fit (`cj_output.txt`, `jwl_output.txt`): eJWLC parameters at the active-products reference and CJ density, `E0` so that the isentrope energy at the transfer pressure equals the JANAF sensible energy of the freeze-out gas at its Abel–Nobel temperature ($e_s = \sum_k Y_k(H_{s,k} - R_k T)$, the blastFoam convention), `omegaT`, the programmed-burn `e0` and the transfer `targetFractions`. Reproduces the previous products-only values (T_trans 825.76 K, omegaT 0.355493); its `E0` is 5.5 kJ/kg below the earlier value, which had left the 0D gas 4.9 K above the products.
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
  - Temperature fix: new optional eJWLC exponent `omegaT` (Section 1.2), keeping the EuDetCode `T_CJ`. `omegaT` was computed with `compute_omegaT.py` and set in all eJWLC cases: 0.355493 (products-only Comp-B/Al fit, 0.371065 with the 2026-10-13 EuDetCode fit: 0D tests, `spherical_alCombustion`, `blastFoam/ChargeCompBAl_1D`), 0.375928/0.375955 (full Comp-B/Al fit with 27% inert Al: `spherical`, `spherical_reacting`, blastEulerFoam tutorials; replaced by the products-only fit on 2026-10-13), 0.469442 (Comp-B into air: `planar_znd`, `planar_znd_phaseTransfer`), 0.353900 (`ChargeCompBAl_1D_AlInert`).
  - The 0D tests now run blastEulerFoam (independent phase energies, so the shared-temperature effect no longer applies): 0.58% maximum pressure deviation, gas formed at 830.6 K, momentum exactly zero, energy conserved to about 5e-9 (since the conservative decode; before it, the dilute-phase reset at 91 µs removed 1.8% / 1.1%).
- **`pressurePhaseChange` in blastEulerFoam moved mass without momentum or energy**: the source kept its momentum and energy while losing mass, and the targets received mass with neither. Now the transferred mass carries the source velocity and total energy (Section 2.1).
- **`tests/0D/phaseTransfer_0D_3phase`** (107.9 bar initial pressure in the blastFoam version, from the mixture-energy closure with an energy-independent Murnaghan phase): in blastEulerFoam the phases have their own energies; the case starts at 100 bar and its pressure check passes (0.58%).
- **Dangling `tmp` references** in `RanzMarshallHeatTransfer.C` and `sphericalHeatTransfer.C` (bound a reference to the temporary returned by `d()`).
- **EuDetCode paths**: the test scripts and the fitter pointed at a non-existent `modules/eudetcode`. They now use `external/eudetcode`, a separate git checkout not tracked by this repository. `tests/1D/verify_planar_znd.py` passes against it (detonation velocity within 0.076%). (2026-10-13) With the EuDetCode GoogleTest binary built, Tier 0 passes (169 tests).
- **0D verifier** (`verify_0d_phase_transfer.py`): energy conservation is checked on the datum-corrected energy, stored energy minus `eOffset` × transferred mass.

- **Double counting of Al** (2026-10-13): the Comp-B/Al cases with a separate Al phase used the full eJWLC fit, which includes the inert Al, and α_Al = 0.3; now the products-only fit and α_Al = 0.182 (Section 6.2).
- **Al phase temperature** (2026-10-13): the granular Al cooled behind the detonation (to the 100 K clip) because its energy equation lacked the compression work and advected kinetic energy into the thermal energy (Section 3.3). It now stays at 288 K under compression and reaches at most 426 K in 40 µs (inert case).
- **Velocity relaxation energy** (2026-10-13): `VRelaxation instant` created or destroyed the dissipated energy and heated trace phases by their own velocity defect; now conservative (Section 3.3). Drag work evaluated at the update-averaged interfacial velocity.
- **EuDetCode CJ state** (2026-10-13): fits regenerated with the current EuDetCode (Section 6.2).
- **Unreacted charge/air interface instability** (2026-10-14): a trace of a dense phase took its own impedance in the HLLC flux, so the charge surface pulled the first air cell back with the explosive's impedance whenever that cell moved outwards; now the trace side uses the mixture impedance (`alphaTrace`, Section 4.6). Tried first: a bracketed search in `equilibrateDiluteVolumes()` for trace explosive stuck in tension at the Murnaghan $p = 0$ limit. It made the interface event later but larger (2.9 kbar at 20 µs), and its trial states could compress the gas past the Abel-Nobel covolume (the inert case failed at 8 µs with maxCo 0.19). Not adopted.
- **gunIB1D comparison** (2026-10-13): with the rebuilt gunIB1D, `verify_gunib1d_benchmarks.py` passes. It finds the binary (`GUNIB1D_DIR`, `build/gunIB1D` or `gunIB1D`) and the cases (`Verification/`, formerly `Validation/`), runs gunIB1D in a removed temporary directory, locates its blastEulerFoam cases next to the script, and first checks Houim 3 against the Houim (2016) reference (compaction speed 0.81%, particle mass to 4e-12). Cross-solver differences: Houim 1 0.000%; Houim 2 shock speed 1.20%, post-shock pressure 0.42%, gas velocity 0.67%; Houim 3 peak pressure 0.59%.

- **Non-conservative `fluidPhaseModel` / `granularPhaseModel` decode** (2026-10-12): dilute phases ($\alpha < 0.01$) had their density, velocity and energy reset with the conserved variables rewritten. This destroyed 1.8% / 1.1% of the energy in the 0D tests and deleted the energy of the products expanding into the air in the spherical cases. Removing the resets exposed two causes of spurious GPa pressures that the resets had masked:
  - the granular buoyancy term built from the individual fluid-phase pressures (Section 3.3);
  - trace fluid phases at fluid/fluid contacts whose decoded density drifts (isostrain without pressure relaxation; fixed by `equilibrateDiluteVolumes()`).
  Both decodes are now conservative. 0D energy is conserved to ~5e-9, and the spherical verifiers (0.64% each since 2026-10-13), the new interface test and the blastEulerFoam tutorials pass. With the gas density floor also removed (2026-10-13), the spherical air shock reaches 10 cm at 17.9 µs with 14 MPa (originally 25.4 µs and 90 MPa). Details, rejected approaches and remaining points: `report/dilute_phase_conservation_plan.md`.

### Open
- **Pressure non-equilibrium in non-dilute mixed fluid cells**: without pressure relaxation, fluid phases with $\alpha > 0.01$ in the same cell undergo the same volumetric strain. At the products/air contact about 50 cells are mixed, with the products at 2–5000 times the gas pressure and the entrained air rarefied at its 5000 K ceiling (2026-10-13). The upstream `instant` pressure relaxation collapses the time step; extending `equilibrateDiluteVolumes()` to all phases (`alphaDilute 1`) equalises the pressures but leaves trace-mass products at 1e5 K. Adding products/gas temperature relaxation fixes this in the inert and reacting cases (Section 6.2); used in `_alCombustion` too since the 2026-10-14 interface fix, and the common temperature in the planar tutorial's mixed cells (about 9000 K) is high. See `report/dilute_phase_conservation_plan.md` (Section 6).
- **`tutorials/blastEulerFoam/reactingParticles`**: after the `d0` fix it runs past 2.6e-7 s, but the time step becomes erratic within the first few hundred steps (down to 1e-15 s) and collapses to 1e-22 s near 2.7e-4 s (of 2.5e-3 s), followed by a floating-point exception in the gas thermo. The original decode behaves the same (collapse at 2.65e-4 s). Not analysed.
- **Verifiers do not check that the solver finished**: `verify_spherical_compBAl*.py` report success from the time directories written before a crash (seen with the ghost-state HLLC); `bin/run_blast_tests.py` checks the log, the scratch runs did not.
- **blastFoam Comp-B/Al tutorials not updated to the new fits**: `tutorials/blastFoam/ChargeCompBAl_1D` has no Al phase but uses the products-only fit (with 1820 kg/m³ reactants and an `e0` entry the eJWLC does not read); `ChargeCompBAl_1D_AlInert` uses an older products fit (see below). The 2D `tutorials/blastEulerFoam/ChargeCompBAl` (8 processes, AMR, 200 µs) was updated but not run.
- **`tutorials/blastFoam/ChargeCompBAl_1D_AlInert`** sets `E0 = -5.484e6` on the products and `eOffset = 5.484e6` in `fvModels`. The two cancel, so the transferred energy is as if `E0 = 0`: the datums are not reconciled for its `hConst` gas. It should use a calibrated `E0` and no `eOffset`.
- **blastFoam mixture closure**: in blastFoam two-/multiphase systems each phase EOS is evaluated with the mixture energy (Allaire-type closure, exact for Mie–Grüneisen phases only); a phase whose pressure is independent of energy (Murnaghan) is then inconsistent. Not relevant to the blastEulerFoam problems.
