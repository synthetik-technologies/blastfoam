# blastfoam

_blastFoam_ is a solver for multi-component compressible flow with application to high-explosive detonation, explosive safety and air blast. 


## Get the code

blastFoam is free. blastFoam is currently in pre-release beta. If you'd like to request a copy of the code and solver before the official public release, while we're in beta, just get in touch: blastfoam@synthetik-technologies.com


## Features

_blastFoam_ currently supports the following features:

- An arbitrary number of phases/EOS's
- JLW equation of state with constant, linear, and "Miller" afterburn models
- Multiple example and tutorial cases
- Automatic mesh refinement (AMR)
- Single and multi-point detonation
- High-order (1st, 2nd, and 4th order in time; 3rd order spatial)
- HLLC, AUSM+, Kurganov, Tadmor flux schemes
- Parallel (MPI)
- Compatible with all standard OpenFOAM mesh generation, pre- and post-processing utilities



## Equations of State

blastFoam includes the following equations of state:

- Jones Wilkens Lee (with afterburn)
- Ideal Gas
- Stiffened Perfect Gas
- Cochran-Chan
- Tait



## Verification and Validation

_blastFoam_ has been validated against known solutions to standard gas dynamics problems, and against data from physical tests.  



## Citation
If you use this code for your work or research, please cite this repository:

```
@software{heylmun_blastfoam:_2019,
	title = {{blastFoam}: An {OpenFOAM} Solver for Compressible Multi-Fluid Flow with Application to High-Explosive Detonation},
	url = {https://github.com/synthetik-technologies/blastfoam},
	publisher = {Synthetik Applied Technologies, {LLC}.},
	author = {Heylmun, Jeffrey and Vonk, Peter and Brewer, Timothy},
	date = {2019-10-22}
}
```
