<!--- python3 -m readme2tex --nocdn --pngtrick --htmlize --output README.md INPUT.md && open README.md.html -->


<p align="center">
  <img src="media/synthetik-logo.png" width="350" title="hover text">
</p>



# blastFoam Version 6.0

blastFoam is a library for single and multiphase compressible flow with application to high-explosive detonation, explosive safety and airblast, as well as general compressible flows. blastFoam is developed by [Synthetik Applied Technologies](https://www.synthetik-technologies.com). This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM<img alt="$\textregistered$" src="svgs/6abda71802c3922eebfcf1b67d5169b2.png" align="middle" width="16.438455000000005pt" height="22.831379999999992pt"/>  and OpenCFD<img alt="$\textregistered$" src="svgs/6abda71802c3922eebfcf1b67d5169b2.png" align="middle" width="16.438455000000005pt" height="22.831379999999992pt"/> trade marks.


## How to use blastFoam

Several validation and tutorial cases are included in the repository, and are documented in the [blastFoam User Guide](BlastFoam_User_Guide.pdf).


### blastFoam workshop | Date: July 14, 2020 | Location: Virtual/Online | Cost: Free

Thank you to everyone who attended!

A recording of the workshop is also available on YouTube: [blastFoam July 14, 2020 Workshop on YouTube](https://youtu.be/OLpRcjR3cW4).


### blastFoam workshop | Date: May 13, 2020 | Location: Virtual/Online | Cost: Free

Thank you to everyone who attended!

The workshop tutorial case and modifications have been added to the repository: [tutorials/blastFoam/building3DWorkshop](tutorials/blastFoam/building3DWorkshop).

A recording of the workshop is also available on YouTube: [blastFoam May 13, 2020 Workshop on YouTube](https://www.youtube.com/watch?v=0YYrkwNWM6U).

The slides from the workshop have been added to the repository: [blastFoam May Workshop Slides](media/Synthetik-Applied-Technologies-blastFoam-Workshop-May-Final-20200513.pdf)


## Installation

Detailed instructions on how to install and use blastFoam are found in the [blastFoam User Guide](blastFoam_User_Guide.pdf). Installation is simple and requires only OpenFOAM-9 and (optionally) gnuplot be installed. Basic installation steps are as follows:



### How to install OpenFOAM for Linux
Compiling OpenFOAM is straight forward, and a more detailed guide to installation can be found [here]{https://openfoam.org/download/source/software-for-compilation}. Once the necessary dependencies have been installed

1. Clone the OpenFOAM-9 repository
```bash
cd $HOME/OpenFOAM
git clone https://github.com/OpenFOAM/OpenFOAM-9.git
```

2. Compile OpenFOAM
```bash
cd OpenFOAM-9
echo "source $HOME/OpenFOAM/OpenFOAM-9/etc/bashrc" >> ~/.bashrc
source etc/bashrc
./Allwmake > log.Allwmake 2>&1
```


### How to install OpenFOAM for Windows 10
An installation video for Windows 10 is available on our YouTube channel: https://youtu.be/vfd610LadSU

<p align="center">
  <img src="media/installation_thumbnail.png" width="450" title="Installing blastFoam in Windows 10">
</p>




### How to Install OpenFOAM for macOS
Compiling OpenFOAM on macOS is relatively straightforward. This [guide and repository](https://github.com/mrklein/openfoam-os-x/wiki/OpenFOAM(R)-git-version-&-Homebrew) provides step-by-step instructions as well as the necessary patch to compile OpenFOAM on macOS.



### How to install blastFoam

1. Install OpenFOAM-9 (if not already installed, see above)
See https://openfoam.org/version/9 for OpenFOAM installation instructions.

2. Create the OpenFOAM directory
```bash
mkdir -p $HOME/OpenFOAM
```

3. Go to the $HOME/OpenFOAM directory
```bash
cd $HOME/OpenFOAM
```

4. Clone the blastFoam repository
```bash
git clone https://github.com/synthetik-technologies/blastfoam.git
```

5. Go to the blastfoam directory
```bash
cd $HOME/OpenFOAM/blastfoam
```

7. Append the etc/bashrc to your .bashrc and/or .zshrc file
```bash
echo "source $HOME/OpenFOAM/blastfoam/etc/bashrc" >> $HOME/.bashrc
# or if using zsh:
echo "source $HOME/OpenFOAM/blastfoam/etc/bashrc" >> $HOME/.zshrc
```


8. Load and set the bash environment to compile blastFoam
```bash
source $HOME/.bashrc
# or if using zsh:
source $HOME/.zshrc
```

9. Compile blastFoam (for parallel use "-j")
```bash
./Allwmake -j
```

10. Test your installation by running the tutorial and validation cases




### Questions
If you find any bugs, please let us know in the issues section of the repository. If you want to get in touch: info@synthetik-technologies.com.




### Citation
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
J. Heylmun, P. Vonk, and T. Brewer, "blastFoam 6.0 User Guide", Synthetik Applied Technologies, LLC., 06-Aug-2020.
```
BiBTex:
```
@misc{heylmun_blastfoamguide_2021,
	title = {{blastFoam version 6.0} {User} {Guide} },
	url = {https://github.com/synthetik-technologies/blastfoam},
	language = {English},
	publisher = {Synthetik Applied Technologies, LLC.},
	author = {Heylmun, Jeffrey and Vonk, Peter and Brewer, Timothy},
	month = aug,
	year = {2022}
}
```






### blastFoam Publications:

## Journals/Conferences

1. T. Brewer, J. Heylmun, and P. Vonk, "Employment of the Open-source Airblast Solver blastFoam to Support the Super Heavy Improvised Explosive Loading Demonstration (SHIELD) Test Program" presented at the ISIEMS, USA, 2019.
2. D. Stephens, P. Vonk, and T. Brewer, "Validation of Open-source Airblast Solver (blastFoam) in an Urban Environment," presented at the MABS 25, Hague, Netherlands, 2018.
3. P. Vonk, "A New OpenFOAM Solver for Compressible Multi-Fluid Flow with Application to High-Explosive Detonation," presented at the OpenFOAM Users Conference, Cologne, Germany, 2016.
4. P. Vonk, T. Brewer, "A New OpenFOAM Solver for Compressible Multi-Fluid Flow with Application to High-Explosive Detonation and Extended Validation," presented at the OpenFOAM Users Conference, USA, 2016.





### blastFoam Version 6.0 Release Notes and Features
blastFoam 6.0 greatly improves the numerics library to include multivariate root finding, minimisation/optimisation, and numerical integration. The equation structure has also been improved and generalised. Lookup tables have also been greatly improved by adding support for mixed order interpolation and more general file readers. 3D lookup tables have also been implemented.

Additional thermodynamic and transport models including power law and  tabulated. Two and multiphase temperature calculations have been expanded to now solve a coupled root finding method that incorperate all phase contributions of energy and density.

Improved and additional refinement methods including polyhedral refinement from foam-extend with the support for 3D, 2D, and 2D axisymmetric refinement with dynamic load balancing.

Improvements for the rotateFields utility now allow for iterative refinement based on the methods provided in dynamicMeshDict.

Burst patches have been added to allow for the breaking of internal faces based on a given criteria.

The setRefinedField utility has been expanded to include the option to set zones and sets based on the given regions.

Time integration now supports restarting of time steps which is useful for FSI cases.


### Previous Releases

### blastFoam Version 5.0 Release Notes and Features
blastFoam 5.0 greatly improves compatibility of blastFoam thermodynamics with that of standard OpenFOAM. This results in the ability to use most standard OpenFOAM functionObjects and fvModels and constraints. This also includes the ability to compile new combinations of thermodynamic models at run-time. Additional thermodynamic models have been added including ePower, ePolynomial, eTabulated, hPower, hPolynomial, and hTabulated. Additional fluid transport models have been added including polynomial, logPolynomial, sutherland, WLF, and tabulated. Additional sold transport models have been added including exponential and polynomial.

Support for fvModels and fvConstraints has been added. This allows for sources that are not typically included in the solver, for example point mass sources.

Improvements for the rotateFields utility included mapping of fields not included in the target directory as well as fixes to the rotation of non scalar fields. 1-D cases can now be directly mapped to 3D cases.

A large selection of numerical methods have also been added including univariate and root finding methods, numerical integration, and minimization.

New methods to initialize non-uniform fields have been added to the setRefinedFields utility such as calculatedDensity, function, and massIntegrate. New cell set types have also been added and include boxMassToCell, cylindericalMassToCell, and sphericalMassToCell.

Because the output of lagrangian particles is not compatible with the builtin paraview reader, the convertLagrangianPositions utility has been added to convert the standard lagrangian format to one that paraview can read, allowing the viewing of parallel largrangian cases, without the need to reconstruct cases.

### blastFoam Version 4.0 Release Notes and Features
blastFoam 4.0 introduces particle solvers. This includes Eulerian multi-fluid methods (blastEulerFoam), Lagrangian (blastParcelFoam), and possible coupling to OpenQBMM for number-density function transport coupling (blastPbeTransportFoam, blastUncoupledVdfTransportFoam, and blastVdfTransportFoam).

Additional equations of state have been added including the Tillotson equation of state for cavitating fluids and the Abel-Nobel equation of state for propellants.

Improvements for activation model including the optional use of delayed detonations and the size of activation points for all activation models, as well as the new programmedIgnition activation model.

Diameter models have been added for use with the Arrhenius rate activation model, as well as the phase models used in the blastEulerFoam solver. These include constant diameter, constant mass, and a quadrature-based method of moments (qbmm) diameter.

New numerical schemes have been added which include 1st, 2nd, and 3rd order MUSCL reconstruction with limiters. Additional options for Strong-stability-preserving Runge-Kutta methods have been added including one, two, and three step RK1-SSP, two, three, and four step RK2-SSP, and three and four step RK3-SSP.

A new functionObjects including tracer particles for the transport of passive particles and the monitoring of conserved quantities.

A simple fluid-structure-interaction solver (blastFSIFoam) has been added, however this is still experimental and not yet stable.

Further improvements have been made to the setRefinedFields utility for more control over setting fields and faster convergence of refinement.


## blastFoam Version 3.0 Release Notes and Features

blastFoam 3.0 now includes thirteen equations of state that allow modeling of diverse materials under extreme conditions, with consideration of phenomenologies such as excitation, dissociation and ionization of nitrogen and oxygen in air at higher energies and temperatures, afterburn, and sympathetic detonation.

We have introduced several different approaches to model detonation within explosive materials which transition from unreacted energetics to detonation products, including pressure-based activation models with multi-step Arrhenius reaction rates, and simple, yet practical models based on empirically derived detonation velocities. Users can also specify instantaneous activation.

blastFoam allows phenomena such as size effect (decrease of the detonation velocity with decreasing charge radius), and detonation front curvature (induced by edge lag of the front as energy is lost to the exterior of the charge) to be accurately captured.  These additions greatly enhance timing accuracy and load characterization, especially for near-contact explosive scenarios.  Options for modeling afterburn (i.e., under-oxygenated explosives continuing to burn after detonation) are also included using the Miller extension, constant, and linear rate models.

blastFoam extends OpenFOAM's base AMR library, and includes the ability to perform 2D and 3D adaptive mesh refinement (AMR). The refinement criteria can be based on density gradient, change across faces (delta), or Lohner's method (2nd derivative of a field) to determine what cell should be refined or unrefined. Additionally, options for mesh unrefinement/relaxation/coarsening have been added, and this is useful for keeping cell counts relatively constant during a calculation while still capturing key features (e.g. shocks) with high accuracy. This allows blastFoam to solve engineering-scale simulations at an affordable computational cost.

blastFoam extends OpenFOAM by adding dynamic load rebalancing for adaptive grids, and now includes a working solution for 2D and experimental support for 3D calculations.  Essentially, at a predetermined timestep interval the domain is rebalanced so that the cell count per CPU is more evenly distributed. This mitigates potential memory issues such as crashing and slow-down related to overloading CPUs that are operating on zones of high refinement.

Turbulence and radiation models have been integrated, allowing blastFoam users to leverage the extensive OpenFOAM libraries and apply them to their simulations, and a new fluid model structure (fluidThermo class), that extends OpenFOAM's standard thermo classes has been added, and provides thermodynamically consistent solutions for more accurate temperature calculations.

New *functionObjects* have been added to improve usability, including the ability to calculate peak overpressure and impulse for each cell in the domain, as well as *blastToVTK*, a utility to view time series mesh surface outputs in *ParaView*.

Additional validation and tutorial cases are also provided to demonstrate and showcase the new functionality and capabilities of blastFoam v3.0.

The code contains multiple utilities to prepare calculations for complex geometries of interest (e.g. engineering-scale; from CAD models), including parallel mesh generation, mesh refinement, advanced post-processing, and import/export functions. Verification and validation studies have been conducted with independent validation (conducted by others) performed on larger-scale problems with complex geometries and published in peer-reviewed journals.  The solver can be run on any modern platform (e.g. laptop, workstation, HPC, AWS, GCP, etc.).

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
- Performance improvements for axisymmetric cases that use wedge boundary conditions (v3.0.1)

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



### Verification and Validation

blastFoam has been validated against known solutions to standard hydrodynamics problems, and against data from physical tests. Validation cases are included with example/tutorial cases as part of the solver source code.




### Disclaimer

This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM<img alt="$\textregistered$" src="svgs/6abda71802c3922eebfcf1b67d5169b2.png" align="middle" width="16.438455000000005pt" height="22.831379999999992pt"/>  and OpenCFD<img alt="$\textregistered$" src="svgs/6abda71802c3922eebfcf1b67d5169b2.png" align="middle" width="16.438455000000005pt" height="22.831379999999992pt"/> trade marks.
