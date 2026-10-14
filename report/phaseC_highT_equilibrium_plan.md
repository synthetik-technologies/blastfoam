# Phase C: High-Temperature Gas Thermodynamics and Equilibrium Chemistry

Plan for step C of `status_review_2026-10-14.md` (WP 1.3 stage 2), extended on 2026-10-14 to include local chemical equilibrium.

## 1. Why

- **Shocked air is too hot.** The gas has no dissociation: the Mach 14.5 air shock test gives T₂ ≈ 10 000 K where equilibrium air gives about 7700 K (`tests/1D/shock_dissociating_air`, fails by design). In the spherical Comp-B/Al cases the shocked air reaches 16 500 K, and the first air cell reaches the 50 000 K decode cap at breakout.
- **The finite-rate chemistry fails at high temperature.** The one-step irreversible reactions (CH₄, CO, H₂ and Al oxidation) are stiff, and the ODE integrators (RKCK45, then seulex) failed with a floating-point exception in the Arrhenius rate in the spherical cases: at breakout (cells at the temperature cap) and at 1.12 ms in the 2 m case (gas at 3360 K). The temporary guard `TreactMax` (no reactions above 6000 K, added 2026-10-14) is not a physical model: it turns the chemistry off where it is fastest.
- At these temperatures the chemical time scales are much shorter than the flow time step, so the gas is close to **local chemical equilibrium**. Equilibrium at high temperature needs the dissociated species (O, N, NO, H, OH, …); with only the major species an equilibrium solve would drive everything to CO₂ and H₂O.

## 2. Approach

A gas cell is treated in one of two ways:

- **Equilibrium** where the gas is hot (above `Teq`, e.g. 2500 K) or the chemistry is fast compared with the time step: the composition is set to chemical equilibrium at the cell's internal energy and density. This conserves mass, elements and energy exactly.
- **Finite rate** elsewhere: kinetics (reversible, so that it tends to the same equilibrium).

## 3. Steps

| # | Step | Exit criterion |
|---|---|---|
| C1 | Species data: converter from the EuDetCode NASA-9 database (`external/eudetcode/Data/GasFlame.json`: three ranges, 200–1000–6000–20 000 K) to OpenFOAM species dictionaries; Al gas species (AL, ALO) from the NASA-7 data (to 6000 K) with constant-c_p extrapolation; condensed Al₂O₃ and C(gr) | Converted c_p, h, s equal the database within round-off; species set N₂, O₂, NO, N, O, H₂, H, OH, H₂O, CO, CO₂, C, CH₄, NH₃, C₂H₆, Al, AlO, Al₂O₃(s), C(gr) |
| C2 | `nasa9BlastThermo`: N-range NASA-9 thermo (constant c_p above the last range), instantiated for the multicomponent/reacting gas phases with Abel–Nobel; the existing `janaf9BlastThermo` has two ranges and is not instantiated for gas phases | 0D test: c_p, h, s of each species and of air against the database to 20 000 K |
| C3 | Equilibrium solver at fixed (e, ρ) (Gibbs minimisation by element potentials, as in CEA; ideal-gas chemical potentials with the Abel–Nobel covolume in the EOS; condensed species as pure phases) | 0D: equilibrium composition and temperature of air and of the Comp-B/Al products against the reference (EuDetCode equilibrium or Cantera C++) to 20 000 K, within 0.5% |
| C4 | Hybrid chemistry in the gas phase: equilibrium above `Teq` (or where τ_chem ≪ Δt), finite rate below; remove `TreactMax` | Elements, mass and energy conserved to round-off in 0D; no ODE failures in the spherical and 2 m cases |
| C5 | Air shock test with equilibrium air | `tests/1D/shock_dissociating_air` passes (jump within 2.5% of equilibrium Rankine–Hugoniot at Mach 14.5) |
| C6 | Finite-rate mechanism below `Teq`: reversible reduced CO/H₂/O₂ kinetics consistent with the equilibrium thermodynamics (plan step B4) | Ignition delay / burnout vs reference mechanism |
| C7 | Re-run spherical cases (0.2 m, 2 m) | Shocked air ≤ ~10 000 K; no temperature-cap hits; documented changes in blast metrics |

## 4. Open choices

- The equilibrium criterion: a fixed `Teq`, or a comparison of the chemical and flow time scales per cell.
- Non-ideality: the equilibrium composition with the Abel–Nobel covolume (affects high-density products near the transfer pressure, b·ρ ≈ 0.05) or ideal-gas chemical potentials.
- Al chemistry: Al vapour and AlO in the equilibrium set, and how the particle burn rate couples to it (plan step D).
- Reference for C3: the EuDetCode equilibrium solver (same data) or Cantera (C++ installed under `/usr/local`; the Python module is not installed).
