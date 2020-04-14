<!--- python3 -m readme2tex --nocdn --pngtrick --htmlize --output README.md INPUT.md && open README.md.html -->


<p align="center">
  <img src="media/synthetik-logo.png" width="350" title="hover text">
</p>



# blastFoam Version 3.0

blastFoam is a solver for multi-phase compressible flow with application to high-explosive detonation, explosive safety and airblast, as well as general compressible flows. blastFoam is developed by [Synthetik Applied Technologies](https://www.synthetik-technologies.com).



## How to use blastFoam

Several validation and tutorial cases are included in the repository, and are documented in the [blastFoam User Guide](blastFoam_User_Guide.pdf).


### blastFoam workshop | Date: May 13, 2020 | Location: Virtual/Online | Cost: Free

An introduction to blastFoam - a free and open-source Computational Fluid Dynamics (CFD) air blast code for modeling high-explosive detonations and blast load generation suitable for blast engineering and protective design.

At the conclusion of the workshop, attendees will be able to competently and independently setup, calculate, visualize and post-process CFD solutions to air blast problems in complex geometries using blastFoam.

The workshop will include:

- An overview of the background theory, physics, equations and implementation in blastFoam
- How to set up and run air blast calculations in simple and complex geometries
- Mesh generation using snappyHexMesh and blockMesh for use with blastFoam
- Setting initial conditions (e.g. JWL EOS coefficients, detonation points, charge shapes, etc.)
- Boundary conditions
- Using probes to get pressure-time histories at discrete locations
- Generating pressure and impulse loads suitable for engineering/protective design
- Pre- and post-processing of results
- Verification and Validation
- Using the blastFoam GUI to help setup and visualize blastFoam cases
- A hands-on guided walk-through of a number of 2D and 3D blastFoam tutorial examples with the Synthetik team.
- All levels of expertise are welcome to attend, and this session will be very practical and hands-on.

Looking forward to connecting with current and new blastFOAMers! Please do reach out with any questions, suggestions, or topics to cover during the workshop and we shall endeavor to work them into the program.

More information and registration here: https://www.eventbrite.com/e/blastfoam-workshop-tickets-100310659884



## Installation

Detailed instructions on how to install and use blastFoam are found in the [blastFoam User Guide](blastFoam_User_Guide.pdf). Instalation is simple and required only OpenFOAM-7 and (optionally) gnuplot be installed. Basic installation steps are as follows:

1. Create the OpenFOAM directory
```bash
mkdir -p $HOME/OpenFOAM
```

2. Go to the $HOME/OpenFOAM directory
```bash
cd $HOME/OpenFOAM
```

3. Clone the blastFoam repository
```bash
git clone https://github.com/synthetik-technologies/blastfoam
```

4. Go to the blastfoam directory
```bash
cd $HOME/OpenFOAM/blastfoam
```

5. Append the etc/bashrc to your .bashrc file
```bash
echo "source $HOME/OpenFOAM/blastfoam/etc/bashrc" >> $HOME/.bashrc
```

6. Load and set the bash environment to compile blastFoam
```bash
source $HOME/.bashrc
```

7. Compile blastFoam (for parallel use "-j")
```bash
./Allwmake
```

8. Test your installation by running the tutorial and validation cases




## Questions and Availability
If you find any bugs, please let us know in the issues section of the repository. If you want to get in touch: info@synthetik-technologies.com. blastFoam is also available on the Texas Advanced Computing Center https://www.tacc.utexas.edu (TACC) as well as several other HPC centers.




## Citation
If you use this code for your work or research, please use this citation:

```
blastFoam: An OpenFOAM Solver for Compressible Multi-Fluid Flow with Application to High-Explosive Detonation. Synthetik Applied Technologies, LLC., 2020.
```
BiBTex:
```
@software{blastfoam,
	title = {{blastFoam}: A Solver for Compressible Multi-Fluid Flow with Application to High-Explosive Detonation},
	url = {https://github.com/synthetik-technologies/blastfoam},
	publisher = {Synthetik Applied Technologies, {LLC}.},
	date = {2020-04-13}
}
```





## User Guide
To cite the [blastFoam User Guide](blastFoam_User_Guide.pdf).:
```
J. Heylmun, P. Vonk, and T. Brewer, "blastFoam 3.0 User Guide", Synthetik Applied Technologies, LLC., 13-Apr-2020.
```
BiBTex:
```
@misc{heylmun_blastfoamguide_2020,
	title = {{blastFoam version 3.0} {User} {Guide} },
	url = {https://github.com/synthetik-technologies/blastfoam},
	language = {English},
	publisher = {Synthetik Applied Technologies, LLC.},
	author = {Heylmun, Jeffrey and Vonk, Peter and Brewer, Timothy},
	month = apr,
	year = {2020}
}
```






## blastFoam Publications:

### Journals/Conferences

1. T. Brewer, J. Heylmun, and P. Vonk, "Employment of the Open-source Airblast Solver blastFoam to Support the Super Heavy Improvised Explosive Loading Demonstration (SHIELD) Test Program" presented at the ISIEMS, USA, 2019.
1. D. Stephens, P. Vonk, and T. Brewer, "Validation of Open-source Airblast Solver (blastFoam) in an Urban Environment," presented at the MABS 25, Hague, Netherlands, 2018.
1. P. Vonk, "A New OpenFOAM Solver for Compressible Multi-Fluid Flow with Application to High-Explosive Detonation," presented at the OpenFOAM Users Conference, Cologne, Germany, 2016.
1. P. Vonk, T. Brewer, "A New OpenFOAM Solver for Compressible Multi-Fluid Flow with Application to High-Explosive Detonation and Extended Validation," presented at the OpenFOAM Users Conference, USA, 2016.






## blastFoam Version 3.0 Release Notes and Features

blastFoam 3.0 now includes thirteen equations of state that allow modeling of diverse materials under extreme conditions, with consideration of phenomenologies such as excitation, dissociation and ionization of nitrogen and oxygen in air at higher energies and temperatures, afterburn, and sympathetic detonation.

We have introduced several different approaches to model detonation within explosive materials which transition from unreacted energetics to detonation products, including pressure-based activation models with multi-step Arrhenius reaction rates, and simple, yet practical models based on empirically derived detonation velocities. Users can also specify instantaneous activation.

blastFoam allows phenomena such as size effect (decrease of the detonation velocity with decreasing charge radius), and detonation front curvature (induced by edge lag of the front as energy is lost to the exterior of the charge) to be accurately captured.  These additions greatly enhance timing accuracy and load characterization, especially for near-contact explosive scenarios.  Options for modeling afterburn (i.e., under-oxygenated explosives continuing to burn after detonation) are also included using the Miller extension, constant, and linear rate models.

blastFoam extends OpenFOAM's base AMR library, and includes the ability to perform 2D and 3D adaptive mesh refinement (AMR). The refinement criteria can be based on density gradient, change across faces (delta), or Lohner's method (2nd derivative of a field) to determine what cell should be refined or unrefined. Additionally, options for mesh unrefinement/relaxation/coarsening have been added, and this is useful for keeping cell counts relatively constant during a calculation while still capturing key features (e.g. shocks) with high accuracy. This allows blastFoam to solve engineering-scale simulations at an affordable computational cost.

blastFoam extends OpenFOAM by adding dynamic load rebalancing for adaptive grids, and now includes a working solution for 2D and experimental support for 3D calculations.  Essentially, at a predetermined timestep interval the domain is rebalanced so that the cell count per CPU is more evenly distributed. This mitigates potential memory issues such as crashing and slow-down related to overloading CPUs that are operating on zones of high refinement.

Turbulence and radiation models have been integrated, allowing blastFoam users to leverage the extensive OpenFOAM libraries and apply them to their simulations, and a new fluid model structure (fluidThermo class), that extends OpenFOAM's standard thermo classes has been added, and provides thermodynamically consistent solutions for more accurate temperature calculations.

New functionObjects have been added to improve usability, including the ability to calculate peak overpressure and impulse for each cell in the domain, as well as *blastToVTK*, a utility to view time series mesh surface outputs in *ParaView*.

Additional validation and tutorial cases are also provided to demonstrate and showcase the new functionality and capabilities of blastFoam v3.0.

The engineering community needs open, verifiable, validated blast and detonation simulation tools.  Currently available tools: 1) are prohibitively expensive to license to run calculations at the scale and number of CPU cores and nodes required to capture key blast phenomena, 2) do not provide access to the underlying code due to concerns over intellectual property rights, 3) use non-universal file formats for pre- and post-processing, and 4) contain export controlled or distribution limited components.  In response to these limitations, Synthetik Applied Technologies leveraged a widely used opensource CFD library as a foundation upon which to develop a new solver suited for high-explosive detonation modeling and simulation, blastFoam.

Synthetik's solver builds upon the most widely utilized opensource CFD platform available today, and is currently deployed at DoD HPC Centers (e.g., AFRL, ARL, ERDC, Navy, ORS).  The code contains multiple utilities to prepare calculations for complex geometries of interest (e.g. engineering-scale; from CAD models), including parallel mesh generation, mesh refinement, advanced post-processing, and import/export functions. Verification and validation studies have been conducted with independent validation (conducted by others) performed on larger-scale problems with complex geometries and published in peer-reviewed journals.  The solver can be run on any modern platform (e.g. laptop, workstation, HPC, AWS, GCP, etc.).

Synthetik is a formal Texas Advanced Computing Center (TACC) Industry Partner, with access to High Performance Computing (HPC) resources on systems such as the new NSF-funded petascale computing system, Frontera, thus allowing Synthetik to develop and test on state-of-the-art HPC systems at scale.

blastFoam currently supports the following features:

- An arbitrary number of phases and EOS's
- Multiple activation and burn models
- Compatiblity with all OpenFOAM's compressible LES and RANS turbulence models
- Extensive verification and validation
- JLW equation of state with constant, linear, and "Miller" afterburn models
- Multiple example and tutorial cases
- Automatic mesh refinement (AMR)
- Blast specific function object for post-processing
- High-order (1st, 2nd, 3rd and 4th order in time; 2nd and 3rd order spatial)
- HLLC, AUSM+, Kurganov, Tadmor flux schemes
- Parallel (MPI)
- Compatible with all of OpenFOAM's standard mesh generation, pre- and post-processing utilities
- Multiple solvers for high-speed reactive flow and deflatration to detonation transition


## Equations of State

blastFoam includes the following equations of state:

- Ideal gas
- Stiffened gas
- Tait
- Van der Waals
- Landau, Stanyukovich, Zeldovich, and Kompaneets (LSZK)
- Jones Wilkens Lee (JWL)
- Cochran-Chan
- Doan-Nickel
- Jones Wilkens Lee C-Form (JWLC)
- Becker Kistiakowsky Wilson (BKW)
- Benedict Webb Rubin (BWR)
- Murnaghan
- Birch Murnaghan (2nd and 3rd order)
- Tabulated


## Activation models

blastFoam includes the following activation models

- None (instantaneous reaction)
- Multi-point linear activation
- Pressure-based
- Arrhenius rate
- Constant rate



## Afterburn models

blastFoam includes the following afterburn models

- None
- Constant
- Linear
- Miller



## Verification and Validation

blastFoam has been validated against known solutions to standard hydrodynamics problems, and against data from physical tests. Validation cases are included with example/tutorial cases as part of the solver source code.




## Disclaimer

This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM<img alt="$\textregistered$" src="svgs/6abda71802c3922eebfcf1b67d5169b2.png" align="middle" width="16.438455000000005pt" height="22.831379999999992pt"/>  and OpenCFD<img alt="$\textregistered$" src="svgs/6abda71802c3922eebfcf1b67d5169b2.png" align="middle" width="16.438455000000005pt" height="22.831379999999992pt"/> trade marks.
