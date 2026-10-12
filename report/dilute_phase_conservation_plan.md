# Conservative Treatment of Dilute Phases and Material Interfaces in blastEulerFoam

Status: **open** (investigation 2026-10-11/12; no change kept in the code).

## 1. Problem

`fluidPhaseModel::decode()` (`4cb49946`) and `granularPhaseModel::decode()` modify the conserved variables of dilute phases:

| Phase model | Action | Threshold |
|---|---|---|
| `fluidPhaseModel` | density blended towards `residualRho`, `alphaRho` rewritten | $\alpha < 10^{-3}$ |
| `fluidPhaseModel` | velocity blended to zero, `alphaRhoU` rewritten | $\alpha < 10^{-2}$ |
| `fluidPhaseModel` | energy reset to $h_e(288.15\,\text{K})$, `alphaRhoE` rewritten | $\alpha < 10^{-2}$ |
| `granularPhaseModel` | blended velocity/energy/$\Theta$; mass removed where $\alpha > \alpha_\text{max}$ | `residualAlpha`, `alphaMax` |

Mass, momentum and energy are not conserved. In `tests/0D/phaseTransfer_0D` and `_3phase`, the energy reset destroys 1.1–1.8% of the total energy in the last step, when the products fall below 1% volume fraction. The same resets act at every material interface, at the detonation front and in every dilute region.

## 2. Why the resets cannot simply be removed

With the decode made fully conservative (conserved fields modified only where $\alpha \le$ `residualAlpha`), the 0D tests conserve energy to about 5e-9. The 1D spherical Comp-B/Al cases then fail. The detonation front itself still travels at the programmed-burn velocity (about 7600 m/s). The verifier's arrival detection ($p > 100$ bar) is triggered by spurious pressure precursors ahead of the front.

Root cause, from the state of `ChargeCompBAl_1D_spherical` at $t = 1\,\mu$s (cells 97–106 of 400, $r \approx 5$ cm, before the detonation arrives):

- The unreacted solid/air interface at the charge edge disintegrates. The unreacted `compBAl` expands to 600–1240 kg/m³ (Murnaghan pressure at its floor of 0).
- The granular Al there reaches $\Theta = 6\times10^5$–$3\times10^7$ m²/s² (granular pressure 0.7–1.8 GPa), with $T$ pinned at the 100 K / 10000 K limits.
- Trace `compBAl` in the air cells decodes to 2000–2555 kg/m³, giving 1.7–11.6 GPa at $\alpha$ down to 1e-5.
- The mixture pressure $p = \sum_k \alpha_k p_k$ carries these pressures into the air.

Mechanism: in trace cells $\alpha_k$ and $\alpha_k\rho_k$ are both tiny and are transported independently, so $\rho_k = \alpha_k\rho_k/\alpha_k$ is unreliable. A stiff EOS (Murnaghan, eJWLC) or granular kinetic theory turns that error into GPa pressures. All cases run with `pressureRelaxationSolver { type none; }`, so nothing holds the stiff/gas interface in mechanical equilibrium. The original resets mask this.

## 3. Approaches tried (all reverted)

| Approach | Result on the spherical cases |
|---|---|
| Fully conservative decode | precursors (+10% apparent velocity, +24% for `alCombustion`); `reacting` crashes |
| Blending of primitives only (conserved fields untouched) | unstable within 2e-7 s: the fluxes are built from primitives inconsistent with the conserved state |
| Conservative pairwise relaxation of dilute phases ($\alpha < 0.01$ or $10^{-4}$) towards the carrier phase (velocity and temperature) | crashes, detonation velocity off; unexplained 2.5e-5 energy leak in `phaseTransfer_0D_3phase` |
| Dilute fluid phases ($\alpha < 0.01$) take the carrier pressure in decode | precursors persist (−50% to +30%), crashes |
| Instant pressure relaxation (`goodBroyden`) | runs, but $\Delta t$ collapses to about 1e-15 s at Courant number 0.2 (an effective wave speed of order 1e11 m/s) |

The patch of the last state (conservative decode plus dilute-pressure equalisation) is not in the repository.

## 4. Plan

Each step has an exit criterion. The final acceptance criterion is: the conserved variables are modified only where $\alpha \le$ `residualAlpha`; the 0D tests conserve energy to < 1e-8; the three spherical verifiers and the blastEulerFoam tutorials pass.

### Step 1: Minimal reproducer and diagnostics
- New case `tests/1D/stationaryInterface`: unreacted Comp-B/Al (eJWLC reactants with Murnaghan, granular Al) next to air at 1 atm, both at rest, no activation. The exact solution is the initial state.
- Diagnostic (debug switch in `fluidPhaseModel` / `granularPhaseModel`): per-step mass, momentum and energy changed by `decode()`, summed per phase and written to the log.
- Use the diagnostic to quantify the energy destroyed by the resets in the spherical cases and tutorials.
- **Exit:** the reproducer shows the spurious pressures with the conservative decode, and the reset losses are known for each regression case.

### Step 2: Interface initialisation
- Check the initial state at the charge edge: phase pressures, $\alpha$ sum, the granular $\Theta$, and whether `setRefinedFields` creates partially filled cells with inconsistent phase states.
- **Exit:** the initial interface cells are in pressure equilibrium with consistent phase states.

### Step 3: Consistent reconstruction of $\alpha$ and $\alpha\rho$
- Review how the blastEuler flux schemes reconstruct the phase fields at faces. Reconstruct $\alpha_k$, $\rho_k$, $\mathbf{u}_k$ and $p_k$ (or $e_k$) and form the face fluxes from them, instead of limiting $\alpha_k$ and $\alpha_k\rho_k$ independently, so that the face density stays within the neighbouring cell densities.
- Requirement: a uniform pressure and velocity across a material interface are preserved exactly (Abgrall's condition).
- **Exit:** the stationary interface stays at rest with pressure perturbations < 1% of ambient for 10 µs, with the conservative decode.

### Step 4: Pressure relaxation
- Find the cause of the $\Delta t$ collapse with `instant` relaxation: which phase and cell gives the large wave speed, and whether the relaxed state is outside the EOS range.
- Make the instant relaxation robust: bounded phase densities, a scalar Newton iteration on the common pressure instead of Broyden, a fall-back to the unrelaxed state when it fails. Alternatively a finite-rate (`ode`) relaxation with a stiff integrator.
- **Exit:** the spherical cases run with pressure relaxation enabled and without a $\Delta t$ collapse.

### Step 5: Granular phase at interfaces
- Replace the non-conservative $\Theta$ and $\alpha_\text{max}$ handling: excess granular energy goes into the internal energy of the phase; at $\alpha_\text{max}$ the grain density is increased (`alphaRho/alphaMax`) rather than mass removed.
- **Exit:** the granular phase conserves mass and energy in the reproducer and in `phaseTransfer_0D_3phase`.

### Step 6: Remove the resets
- Make both `decode()` functions conservative (the 2026-10-11 version: conserved fields modified only where $\alpha \le$ `residualAlpha`; velocity cap and energy clamp on primitives only).
- Run the full regression (Tier 1 and 2, blastEulerFoam tutorials).
- **Exit:** the final acceptance criterion above.

## 5. Related open items
- `phaseTransfer_0D` WP 1.1 criterion $\Delta p/p < 0.01\%$ (currently 0.58–0.59%).
- Possible double counting of Al in the full Comp-B/Al cases (`UPDATES.md`, Section 7).
