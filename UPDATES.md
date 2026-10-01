# Library Updates (Compared to `origin/develop` @ 71b97fd)

This document provides a summary of updates, new features, thermodynamic models, numerical improvements, and bug fixes introduced in the library relative to `origin/develop` (commit `71b97fd7c9195fdb342ed4230d258dd1dde80b55`).

---

## 1. New Thermodynamic Models & Equation of States

### 1.1 Multicomponent & Reacting Phase Fluid Thermodynamics
- **`multicomponentPhaseFluidBlastThermo`**:
  - Implemented multicomponent phase fluid thermodynamics for multiphase systems (`phaseFluidBlastThermo`), enabling individual phases in a multiphase system to consist of multiple chemical species with species transport and mixture properties.
- **`reactingPhaseFluidBlastThermo`**:
  - Added reacting multicomponent phase fluid thermodynamics, pairing phase multicomponent thermodynamics with `blastCombustionModel` to enable per-phase chemical reactions and combustion modeling in multiphase flow systems.
  - Automatically writes and exports phase-specific heat release rate `Qdot` ($[\text{W/m}^3]$).

### 1.2 `eJWLC` Equation of State
- **Energy-based Mie-Grüneisen JWL with 3-Term Isentrope (`eJWLC`)**:
  - Added the `eJWLC` equation of state implementing the full 3-term ($A, B, C, R_1, R_2, \omega$) isentrope pressure formulation:
    $$p(\rho, e) = A \exp(-R_1 V) + B \exp(-R_2 V) + C V^{-(\omega+1)} + \omega \rho \left(e - e_s(V)\right)$$
    where $V = \rho_0 / \rho$ and $e_s(V)$ is the cold-curve isentropic energy.
  - Implemented an analytic CJ-isentrope reference temperature:
    $$T_s(V) = T_\text{CJ} \left(\frac{V_\text{CJ}}{V}\right)^\omega$$
  - Provides proper EOS energy corrections to recover correct isentropic and off-isentropic temperatures when coupled with `eConst` caloric models.
  - Registered with thermo model selector macros in `forBlastGases.H` and `forDetBlastGases.H`.

---

## 2. Finite Volume Models (`fvModels`) & Phase Change

### 2.1 Pressure-Driven Phase Change Model (`pressurePhaseChange`)
- Added a new `fvModel` for pressure-induced mass transfer from a source phase to one or more multi-component target phases (`pressurePhaseChange`, replacing the legacy `temperaturePhaseChange` model):
  - Driven by a selected pressure interval $[p_\text{lower}, p_\text{upper}]$ (or `pInterval (pMin pMax)`) rather than temperature thresholds, overcoming the well-known inaccuracies of empirical JWL temperature/caloric formulations.
  - Implements a smooth hyperbolic tangent transition factor $\Gamma(p) = \frac{1}{4}\left(1 + \tanh\frac{p_\text{upper} - p}{p_\text{range}}\right)\left(1 + \tanh\frac{p - p_\text{lower}}{p_\text{range}}\right)$ with relaxation timescale $\tau$.
  - Supports target phase species distribution specified via mass fractions or mole fractions.
  - Dynamically calculates the inter-phase internal energy offset between the source EOS and target mixture EOS to strictly conserve total energy to machine precision.
  - Integrates with OpenFOAM's `fvModel` framework and adds conservative source terms to phase continuity equations ($\alpha_i, \rho_i, \alpha\rho_i$) and the total energy equation.

---

## 3. Solver & Compressible System Enhancements

### 3.1 `fvModels` Integration in Compressible Systems
- **`compressibleBlastSystem`**: Added support for energy source terms evaluated from active `fvModels` (`models().source(rhoEff(), e_)`) in `addSources()`.
- **`multiphaseCompressibleSystem`**: Added support for phase mass source terms evaluated from active `fvModels` (`models().source(alphas_[phasei], rhos_[phasei])`) in `solveMass()`.

### 3.2 Two-Phase System Stability & Reconstruction
- **`twoPhaseCompressibleSystem`**:
  - Revamped `decode()` routine with a microscopic volume fraction cutoff ($10^{-15}$) to prevent unphysical EOS evaluation in empty cells.
  - Added density clamping for tiny phase volume fractions ($\alpha < 0.01$) to eliminate localized pressure/density spikes.
  - Improved robustness of phase density and volume fraction reconstruction.

---

## 4. Thermodynamics, Chemistry & Numerics Refinements

### 4.1 Mass Fraction Handling & Mixture Robustness
- **`basicSpecieBlastMixture`**:
  - Overhauled `correctMassFractions()` to gracefully handle empty / zero-density cells where $\sum Y_i \approx 0$, setting the cell composition to the designated `defaultSpecie` (or species 0) without triggering fatal termination.
  - Supported `defaultSpecie` / `inertSpecie` lookup when phases are named.

### 4.2 Chemistry Solvers & Combustion Models
- Generalized `blastChemistryModel` and chemistry solvers (`blastChemistrySolver`, `noBlastChemistrySolver`, `odeBlastChemistrySolver`, `EulerImplicitBlastChemistrySolver`) to support `multicomponentBlastThermo` and `phaseFluidBlastThermo`.
- Instantiated chemistry solvers for `phaseFluidBlastThermo` across standard thermodynamic combinations.
- Updated `reactingFluidBlastThermo` and `reactingSolidBlastThermo` to maintain consistent `Qdot_` heat release rate fields and initialized mass fractions.

### 4.3 Equation of State & Caloric Model Updates
- **`JWL`**: Added configurable maximum pressure cutoff `pMax` (default $10^{15}\,\text{Pa}$) across all constructor overloads and member operators to prevent runaway pressures during extreme compression states.
- **`Murnaghan`**: Added `pMax_` support across constructors and operators.
- **`perfectGas`**: Added maximum pressure bound ($10^{15}\,\text{Pa}$) when pressure limiting is active.
- **`janafBlastThermo`**: Polynomial evaluations for $C_p$, $H_a$, $C_{pt}$, $H_{at}$, $S$, $G_\text{std}$, and $dC_p/dT$ now strictly use temperature-clamped bounds (`limit(T, T)`) to avoid divergence outside valid polynomial ranges.
- **`twoPhaseFluidBlastThermo`**:
  - Clamped denominator `XiSum` ($\ge 10^{-6}$) in mixture pressure calculations.
  - Protected sound speed evaluation against negative roots in degenerate cells.
  - Added diagnostics for NaN and extreme pressure states.

### 4.4 Flux Scheme Selectors
- **`fluxScheme`**: Added overloaded `NewSingle` and `NewMulti` factory methods accepting an explicit `fluxSchemeType` parameter, enabling runtime specification and override of flux schemes by name.

### 4.5 Build System Updates
- **`src/finiteVolume`**: Updated `Make/files` and `Make/options` to compile `fvModels/pressurePhaseChange` and link thermophysical model libraries (`physicalProperties`, `fluidThermophysicalModels`, `specie`, `blastNumerics`, and `thermodynamicModels`).
- **`src/thermodynamicModels`**: Updated `Make/files` to compile `multicomponentPhaseFluidBlastThermos.C` and `reactingPhaseFluidBlastThermos.C`.

---

## 5. Tutorials & Verification Cases

Added and updated cases showcasing multiphase multicomponent flows, reacting flows, and phase change:
- **`tutorials/blastFoam/pressurePhaseChangeTest`**: Demonstrates the `pressurePhaseChange` `fvModel` transitioning condensed phase material to gas products based on pressure intervals.
- **`tutorials/blastFoam/twoChargeDetonation`**: Multi-charge explosive detonation interaction simulation.
- **`tutorials/blastFoam/ChargeCompBAl` / `ChargeCompBAl_1D` / `ChargeCompBAl_1D_AlInert`**: 1D and multi-dimensional aluminized Composition B explosive detonation and afterburning test cases.
- **`tutorials/blastFoam/multicomponent_shock_bubble`**: Shock interaction with a multicomponent gas bubble.
- **`tutorials/blastFoam/multiphase_reacting_bubble` / `multiphase_reacting_bubble_1d`**: Multiphase reacting bubble ignition and expansion.
- **`tutorials/blastFoam/axisymmetricChargeMod` / `axisymmetricChargeReacting`**: Axisymmetric high-explosive detonation with reacting multicomponent products.
- **`tutorials/blastReactingFoam/obliqueShock`**: Updated parallel run scripts and grid configuration.
