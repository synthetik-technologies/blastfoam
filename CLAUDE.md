# blastFoam Project Context

## Project Overview
`blastFoam` is an open-source library and suite of solvers for single and multiphase compressible flow, designed for high-explosive detonation, explosive safety and airblast. Developed by [Synthetik Applied Technologies](https://www.synthetik-technologies.com), it extends OpenFOAM with numerical flux schemes, equations of state (EOS) and adaptive mesh refinement (AMR).

This checkout is a fork (`origin` = `github.com:kenoba/blastfoam`, branch `develop`). It extends blastFoam for aluminised-explosive detonation, phase transfer and afterburning modelling. Upstream `synthetik-technologies/blastfoam` is merged in periodically; the fork's local work starts after `71b97fd7`.

The aluminised-charge problems target **blastEulerFoam** (multi-fluid: per-phase mass, momentum and total energy, with finite-rate velocity/temperature/pressure relaxation). blastFoam's two-/multiphase systems use a single mixture energy and velocity.

### Key Technologies
- **OpenFOAM 12** (`BLAST_OF_VERSION=12`; the top-level `Allwmake` refuses other versions).
- **C++** with OpenFOAM `wmake`; **MPI** for parallel runs; **Python 3 + numpy** for test verification and report tooling.

## Project Structure
- `src/`: libraries, built by `src/Allwmake` in this order: `numerics`, `finiteVolume`, `diameterModels`, `thermodynamicModels`, `radiationModels`, `MomentumTransportModels`, `compressibleSystem`, `polyTopoChange`, `functionObjects`, `regionModels`, `errorEstimators`, `fvMotionSolvers`, `fvMeshTopoChangers`, `fvMeshDistributors`, `atmosphereModels`, `solidModels`, `levelSet`, `FieldSetTypes`, then `regionSolvers/Allwmake` and `solids4Foam/Allwmake`.
    - `numerics` (`libblastNumerics`): generic maths: root finding, minimisation, integrators, lookup tables, least squares, `Function3/4`.
    - `finiteVolume` (`libblastFiniteVolume`): flux schemes (HLLC, AUSM+, Kurganov, Tadmor, … in `fluxSchemes/`), reconstruction, time integrators, cell sources.
    - `thermodynamicModels` (`libblastThermodynamics`): EOS (`specie/equationOfStates/equationOfStates/`: JWL, eJWLC, AbelNobel, Murnaghan, BKW, …), blast thermo classes (`derived/`), activation/programmed burn, afterburn, chemistry and combustion, cavitation/phase change, and `fvModels/pressurePhaseChange`.
    - `compressibleSystem`: flow systems. Base class at `compressibleSystem/compressibleSystem/compressibleSystem.H`, plus `compressibleBlastSystem`, single/two/multiphase(-interface) variants and `coupledSystems`.
    - AMR / load balancing: `fvMeshTopoChangers` (`fvMeshRefiner`), `errorEstimators`, `polyTopoChange`, `fvMeshDistributors`.
    - `functionObjects`: includes `temporalMax` and `temporalIntegral` (peak and impulse tracking).
- `applications/`: built by `applications/Allwmake` in the order `utilities`, then `modules`, then `solvers`.
    - `solvers/`: `blastFoam`, `blastEulerFoam`, `blastParcelFoam`, `blastReactingFoam`, `blastRegionFoam`, `blastXiFoam`, `blastFESolid`, `blastOptimization`, `solidFoam`. Each is a standalone `main()`.
    - `modules/`: solver-module libraries `libblast{Euler,Parcel,Reacting,Solid,Xi}`. Most Euler–Euler multiphase physics lives in `modules/blastEuler` (`phaseSystems`, `interfacialModels` (drag, heat transfer), `kineticTheoryModels`). `modules/blast` has no `Make/` and is not built.
    - `utilities/`: e.g. `setRefinedFields`, `addInternalPatch`, `addEmptyPatch`, `blastToVTK`, `calculateImpulse`, `initializeAtmosphere`, `mergeProbes`.
    - `test/Test-*`: unit-test applications; not built by `Allwmake`.
- `tests/`: the fork's verification cases (`0D/`, `1D/`, `1D_spherical/`), each with a `verify_*.py` script.
- `tutorials/`: example cases. `validation/`: upstream validation cases.
- `bin/run_blast_tests.py`: regression test runner. `utilities/` (top level): `thermo_fitter/` (eJWLC / Abel-Nobel fitter `fit_ejwlc_nobleabel.py`; `compute_omegaT.py` computes the eJWLC `omegaT` for a case; `products_only_fit.py` the Comp-B/Al entries from the EuDetCode output in `eudet_compBAl/`) `mechanism/` (`make_gas_mechanism.py`: reversible GRI-3.0 H2/O2/CO/C1 subset for the gas phase below `Teq`; `cantera_cv_ignition.cpp`: Cantera reference reactor) and `generate_report_figures.py` (report figures from test results under `$BLAST_RESULTS_DIR`, default the repository; run the tests in a copy, e.g. `BLASTFOAM_ROOT=<copy> bin/run_blast_tests.py --tier all`, then point `BLAST_RESULTS_DIR` at it).
- `report/`: `development_plan.md` (staged plan with current status), `stage*_*.md` (per-stage plans), `dilute_phase_conservation_plan.md` (conservative decode and dilute phases), `TBX-blastFoam.tex` (technical report).
- `etc/`: `bashrc` and code templates. `Makefile`, `Dockerfile`, `debian/`: packaging.

## Building and Running

### Environment Setup
```bash
source /opt/OpenFOAM/OpenFOAM-12/etc/bashrc
source /opt/OpenFOAM/blastfoam/etc/bashrc
```
`etc/bashrc` sets `BLAST_DIR` and `BLAST_MODULES`. `BLAST_APPBIN` and `BLAST_LIBBIN` default to `$FOAM_USER_APPBIN` and `$FOAM_USER_LIBBIN`; they are `$BLAST_DIR/bin` and `$BLAST_DIR/lib` only when `MAKE=True` or `/opt/blastfoam` exists.

### Compilation
```bash
./Allwmake -j        # everything
wmake                # inside a single library/application directory
```
- **Library dependency rule:** `libblastThermodynamics` links against `libblastFiniteVolume` and is built after it. Code in `src/finiteVolume` must not use classes from `src/thermodynamicModels`, even though its `Make/options` has the thermodynamicModels include path: it compiles, but the link fails. Put such code in `src/thermodynamicModels` or a later library (see commit `c3294dc1`).
- **Leftover build folders:** after sources are removed or moved, git-ignored `Make/linux64*` folders can remain. `wmake -all` then stops with "file 'Make/files' does not exist". Delete such directories if git tracks no files in them.
- Check the build log for `error`/`undefined` lines. `Allwmake` itself can exit 0 even when sub-builds failed.

### Running a Case
1. Mesh: `blockMesh` (or `snappyHexMesh`).
2. Initial fields: `setRefinedFields` or `setFields` (the 0D tests use `python3 ../setup_case.py`).
3. Solve: `blastFoam` / `blastEulerFoam` / … (serial) or `mpirun -np <N> <solver> -parallel`.

### Tests
```bash
python3 bin/run_blast_tests.py --tier 0|1|2|3|all [--case <relpath>] [--report blast_test_report.md]
```
- The runner runs `Allclean` and `Allrun` for each case and checks the solver log for FATAL, SIGFPE, NaN or a missing `End`. It then runs the matching `verify_*.py` and overwrites `blast_test_report.md`.
- Paths are hard-coded to `/opt/OpenFOAM/blastfoam`.
- `Allrun` scripts delete and regenerate time directories and fields. To try a case without touching the repository, copy it (with its `verify_*.py`, and for `tests/0D` also `setup_case.py`) to a scratch directory. The verifiers locate cases relative to their own directory.
- External dependencies:
    - Tier 0 and the planar ZND verifiers use the EuDetCode build in `external/eudetcode`, a separate git checkout not tracked by this repository (Tier 0 needs its GoogleTest binary `build/tests/eudetcode_tests`).
    - `tests/1D/verify_gunib1d_benchmarks.py` uses a local gunIB1D build.
- The 0D transfer tests (`tests/0D/phaseTransfer_0D`, `phaseTransfer_0D_3phase`) run blastEulerFoam; the 3-phase variant uses `setup_case.py --three-phase` and `verify_0d_phase_transfer.py --three-phase` (both in Tier 1).

## Development Conventions
- **Coding style:** OpenFOAM C++ conventions (PascalCase types, camelCase functions and variables, 4-space indentation, OpenFOAM file headers).
- **Verification:** new physics comes with a case under `tests/` and a `verify_*.py` script, registered in `bin/run_blast_tests.py`. Demonstration cases go in `tutorials/`.
- **Commits:** fork work uses conventional commits (`feat(scope):`, `fix(scope):`, `docs(scope):`). Upstream uses `DEV:`, `BUG:`, `ENH:` and `PRJ:` prefixes.
- **README:** `README.md` is generated from `INPUT.md` (readme2tex; see the comment at the top of `INPUT.md`). Edit `INPUT.md`.
- **Conservation in blastEulerFoam:** interphase transfer must move mass, momentum and energy together. The phase `decode()` functions may modify the conserved variables (`alphaRho`, `alphaRhoU`, `alphaRhoE`) only where a phase is residual; bounds elsewhere act on the primitives only. Do not blend primitives without changing the conserved fields (the fluxes are built from the primitives). Dilute fluid phases are kept near pressure equilibrium by `phaseSystem::equilibrateDiluteVolumes()` (`alphaDilute`), which exchanges only volume and the pΔα work; with `alphaDilute 1` plus a stiff `constantCoefficient` products/gas heat transfer, mixed cells reach pressure and temperature equilibrium (used in the Comp-B/Al cases). Couple phases through the mixture pressure, not through individual phase pressures, which are undefined where a phase is absent. `residualRho` is the reference density of residual cells, not a density floor. Granular packing must be limited by the frictional stress: use a small `alphaDeltaMin` (≤ 1e-3; 0.05 caps the stress at ~10 kPa), since the `alphaMax` bound in `decode()` then stores the excess mass as overdense grains. Trace fluid phases at a face take the mixture impedance in HLLC with `alphaTrace` (`phaseProperties`, 1e-3 in the Comp-B/Al cases; without it a trace of explosive at the charge/air interface acts with the explosive's impedance and the interface goes unstable). `tests/1D_spherical/ChargeCompBAl_1D_spherical_interface` checks this; background in `report/dilute_phase_conservation_plan.md`.
- **Comp-B/Al cases with a separate Al phase** use the EuDetCode products-only eJWLC fit (`--jwl-products-only`, ρ0 = 1624.45) with α_compBAl = 0.818, α_Al = 0.182 (27% inert Al by mass at 1820 kg/m³); the full fit (ρ0 = 1820) already contains the Al. Al burning (`reactingDiameterModel`) is limited to gas-phase oxidiser via `oxidiserPhase`/`oxidiserSpecies`/`YoxidiserRef`. The fit values (current EuDetCode: D = 7154 m/s) and the products `E0`/`omegaT`/initiation `e0` come from `utilities/thermo_fitter/products_only_fit.py`; `E0` uses the blastFoam gas sensible energy $\sum Y_k(h_{s,k} - R_k T)$. The gas is multi-species JANAF / Abel-Nobel air; `janafBlastThermo` extrapolates with constant Cp above the 5000 K fit range (gas `THigh 50000`), without dissociation, so strongly shocked air is too hot. A `nasa9` thermo (three ranges to 20 000 K, species file from `utilities/thermo_fitter/nasa9_species.py`) exists for the gas phases; the gas chemistry sets cells above `Teq` to local equilibrium at fixed (e, ρ) (`equilibrium { Teq 2500; }` in `chemistryProperties`, solver `src/thermodynamicModels/chemistryModel/equilibrium/blastChemicalEquilibrium`; frozen where it fails; finite rate below `Teq`; `report/phaseC_highT_equilibrium_plan.md`). The phase energies are sensible energies: energy checks with gas chemistry must add the gas formation enthalpy and subtract that of the transferred products (as in `tests/1D_spherical/verify_spherical_compBAl_reacting.py`).
- **Keep docs in sync:** after code changes, update `UPDATES.md` (changelog relative to `71b97fd7`) and the Status columns in `report/development_plan.md`.

## Key Files
- `etc/bashrc`: environment variables (`BLAST_DIR`, `BLAST_APPBIN`, `BLAST_LIBBIN`, `BLAST_OF_VERSION`).
- `src/compressibleSystem/compressibleSystem/compressibleSystem.H`: base class for the compressible flow systems.
- `applications/solvers/blastFoam/blastFoam.C`: main solver entry point.
- `src/thermodynamicModels/fvModels/pressurePhaseChange/`: pressure-driven condensed → gas phase transfer.
- `UPDATES.md`: changes made in this fork, with known issues.
- `report/development_plan.md`: stage plan and current status.
- `BlastFoam_User_Guide.pdf`: upstream user and theory guide.
