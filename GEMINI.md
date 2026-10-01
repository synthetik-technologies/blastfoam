# blastFoam Project Context

## Project Overview
`blastFoam` is an open-source library and suite of solvers for single and multiphase compressible flow, specifically designed for applications involving high-explosive detonation, explosive safety, and airblast. Developed by [Synthetik Applied Technologies](https://www.synthetik-technologies.com), it extends the OpenFOAM toolbox with advanced numerical schemes, equations of state (EOS), and adaptive mesh refinement (AMR) capabilities.

### Key Technologies
- **OpenFOAM:** Built on OpenFOAM (specifically version 12).
- **C++:** Core implementation language.
- **MPI:** For parallel execution.
- **Numerical Methods:** High-order schemes (HLLC, AUSM+, Kurganov-Tadmor) and multidimensional AMR.

## Project Structure
- `src/`: Core library code divided into modules:
    - `compressibleSystem`: Main physical system definitions.
    - `thermodynamicModels`: Implementation of various EOS (Ideal Gas, JWL, etc.).
    - `numerics`: Riemann solvers and flux schemes.
    - `fvMeshTopoChangers`: AMR and mesh refinement logic.
- `applications/`:
    - `solvers/`: Primary applications like `blastFoam`, `blastEulerFoam`, `blastParcelFoam`.
    - `utilities/`: Custom pre- and post-processing tools (e.g., `setRefinedFields`, `addInternalPatch`, `blastToVTK`).
- `tutorials/`: Example cases demonstrating various features.
- `validation/`: Cases for verifying solver accuracy against analytical or experimental data.
- `etc/`: Configuration files and environment setup scripts.

## Building and Running

### Environment Setup
Before building or running, you must source the OpenFOAM-12 environment and the blastFoam-specific `bashrc`:
```bash
# Example (adjust paths as necessary)
source /opt/OpenFOAM/OpenFOAM-12/etc/bashrc
source /opt/OpenFOAM/blastfoam/etc/bashrc
```

### Compilation
The project uses standard OpenFOAM `wmake` scripts. To compile the entire project:
```bash
./Allwmake -j
```
Individual solvers or libraries can be compiled by running `wmake` within their respective directories.

### Running a Case
Cases are structured as standard OpenFOAM cases. A typical workflow (as seen in `tutorials/`) involves:
1.  Generating mesh: `blockMesh` or `snappyHexMesh`.
2.  Pre-processing: `setRefinedFields`.
3.  Execution: `blastFoam` (serial) or `mpirun -np <N> blastFoam -parallel` (parallel).

## Development Conventions
- **Coding Style:** Adheres to OpenFOAM's C++ coding standards (e.g., PascalCase for types, camelCase for functions/variables, specific indentation).
- **Scripts:** Use `Allrun` and `Allclean` scripts in case directories for automated testing and cleanup.
- **Verification:** New features should be accompanied by a validation case in the `validation/` directory.

## Key Files
- `etc/bashrc`: Defines environment variables like `$BLAST_DIR`, `$BLAST_APPBIN`, and `$BLAST_LIBBIN`.
- `src/compressibleSystem/compressibleSystem.H`: Central class for the compressible flow system.
- `applications/solvers/blastFoam/blastFoam.C`: Main entry point for the flagship solver.
- `BlastFoam_User_Guide.pdf`: Comprehensive documentation on usage and theory.
