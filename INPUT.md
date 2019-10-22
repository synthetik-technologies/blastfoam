## Disclaimer

This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM$\textregistered$  and OpenCFD$\textregistered$ trade marks.




## What is blastFoam?

blastFoam is a solver for multi-component compressible flow with application to high-explosive detonation, explosive safety and airblast. The blastFoam solver uses OpenFOAM technology, and is in no way approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM$\textregistered$  and OpenCFD$\textregistered$ trade marks (see Disclaimer above).




## Get the code

blastFoam is free and opensource. blastFoam is currently in pre-release beta, and is being made available primarily to academic, research and government. If you'd like to request a copy of the code and solver before the official public release, while we're in beta, just get in touch: info@synthetik-technologies.com




## Citation

If you use this code for your work or research, please use this citation:

```
@software{heylmun_blastfoam:_2019,
	title = {{blastFoam}: A Solver for Compressible Multi-Fluid Flow with Application to High-Explosive Detonation},
	url = {https://github.com/synthetik-technologies/blastfoam},
	publisher = {Synthetik Applied Technologies, {LLC}.},
	author = {Heylmun, Jeffrey and Vonk, Peter and Brewer, Timothy},
	date = {2019-10-22}
}
```





## Features

blastFoam currently supports the following features:

- An arbitrary number of phases/EOS's
- JLW equation of state with constant, linear, and "Miller" afterburn models
- Multiple example and tutorial cases
- Automatic mesh refinement (AMR)
- Single and multi-point detonation
- High-order (1st, 2nd, 3rd and 4th order in time; 2nd and 3rd order spatial)
- HLLC, AUSM+, Kurganov, Tadmor flux schemes
- Parallel (MPI)
- Compatible with all of OpenFOAM's standard mesh generation, pre- and post-processing utilities




## Equations of State

blastFoam includes the following equations of state:

- Jones Wilkens Lee (JWL) (with afterburn)
- Ideal Gas
- Stiffened Gas
- Cochran-Chan
- Tait




## Verification and Validation

blastFoam has been validated against known solutions to standard gas dynamics problems, and against data from physical tests. Validation cases are included with example/tutorial cases as part of the solver source code. 





### Validation/Example: Internal Detonation 

Reference:
```
Joachim, Charles E., Gordon W. McMahon, Christo V. Lunderman, and Sharon B. Garner. 1999. “Airblast Effects Research: Small-Scale Experiments and Calculations.” DTIC Document.
```


Validation against experimental and simulated (CTH) data as reported in Joachim et. al.; see the paper for an explanation of scaling.

![Joachim case setup](media/joachimCaseSetup.png)

![Joachim validation case (scaled)](media/pressureScaledTimePawm.gif)



### Validation/Example: Blast Loads Behind Vertical Walls

Reference:
```
M. E. Beyer, “Blast Loads Behind Vertical Walls,” Naval Civil Engineering Laboratory, Port Hueneme, CA, USA, AD-A181 274, 1986.
```

![Blast Wall Case](media/blastWallwm.gif)





### Validation/Example: Multi-Component Flow Verification

Reference:
```
Zheng, H. W., C. Shu, Y. T. Chew, and N. Qin. “A Solution Adaptive Simulation of Compressible Multi-Fluid Flows with General Equation of State.” International Journal for Numerical Methods in Fluids 67, no. 5 (2011): 616–637.
```

Verification and validation results as compared to those published by Zheng et. al.; HLLC flux shown.


![Verification plots compared with Zheng](media/zheng-blastfoam-validation.png)
















## Governing Equations

The evolution of a two phase, compressible, and inviscid mixture can be defined by a set of coupled evolution equations for mass, momentum, and energy. 

$$
    \partial_t \mathbf{U} + \nabla \cdot \mathbf{F} = \mathbf{S}
$$

$$
    \mathbf{U} = 
        \left( \begin{array}{c}
           \alpha_1 \\
           \alpha_1 \rho_1 \\
           \alpha_2 \rho2 \\
           \rho \mathbf{u} \\
           \rho E
        \end{array} \right)
    \mathbf{F} = 
        \left( \begin{array}{c}
           \alpha_1 \mathbf{u} \\
           \alpha_1 \rho_1 \mathbf{u} \\
           \alpha_2 \rho_2 \mathbf{u} \\
           \rho \mathbf{u} \otimes \mathbf{u} + p \mathbf{I}\\
           (\rho E + p) \mathbf{u}
        \end{array} \right)
    \mathbf{S} = 
        \left( \begin{array}{c}
           \alpha_1 \nabla \cdot \mathbf{u} \\
           0 \\
           0 \\
           0 \\
           0
        \end{array} \right)
$$

where $\rho$ is the mixture density, $\mathbf{u}$ the mixture velocity, $E$ the total energy, $p$ the pressure, and $rho_i$ and $\alpha_i$ are the density and volume fraction of each phase. 





$$
    \alpha_2 = 1 - \alpha_1
$$


$$
    \rho = \sum_i \alpha_i \rho_i
$$


$$
    \rho E = \rho e + \frac{1}{2}\rho |\mathbf{u}|^2
$$

$$
    \rho e = \sum_i \alpha_i \rho_i e_i
$$

The pressure will be defined using a specified equation of state where the mixture internal energy, densities, and volume fraction are used to calculate the total pressure. The equations of state will be used in the Mie-Gruneisen form.

$$
    p_i(\rho_i, e_i, \rho) = (\Gamma(\rho_i) - 1) \rho_i e_i - \Pi(\rho_i)
$$







## Equations of State 


The mixture pressure is defined using the Mie-Gruneisen from using
$$
    p = \frac{\rho e}{\sum_i \alpha_i \xi_i} - \frac{\sum_i \alpha_i \xi_i \Pi_i}{\sum_i \alpha_i \xi_i}
$$ 

where 

$$
    \xi_i(\rho_i) = \frac{1}{\Gamma_i - 1}
$$

and $\Pi_i$ is dependent on the equation of state.

The speed of sound within a give phase is give by

$$
    c_i = \sqrt{\frac{\sum_i y_i \xi_i c_i^2}{\sum_i \xi_i}}
$$

with

$$
    y_i = \frac{\alpha_i \rho_i}{\rho}
$$



$$
    c_i^2 = \frac{h_i - \delta_i}{\xi_i}
$$



$$
    h_i = \frac{\Gamma_i p + \Pi_i}{(\Gamma_i - 1)\rho_i}
$$

and $\delta_i$ is again dependent on the equation of state.












### Generalized van der Waals gas

For a gas described by the generalized van der Waals equation of state, the pressure is defined as
$$
    p_i = \frac{\gamma_i - 1}{1 - b_i \rho_i}(\rho_i e_i + a_i \rho_i^2) - (a_i \rho_i^2 + c_i) 
$$

where $a_i$, $b_i$, $c_i$, and $\gamma_i$ are material parameters. Writing the above equation in M.G. form, we obtain

$$
    \Gamma_i = \frac{\gamma_i - 1}{1 - b_i \rho_i} + 1 
$$

$$
    \Pi_i = 
        \left[1 - \frac{\gamma_i - 1}{1 - b_i \rho_i}\right] a_i \rho_i^2 
      + \left[\frac{\gamma_i - 1}{1 - b_i \rho_i} + 1\right] c_i 
$$

and

$$
    \delta_i = 
      - b_i \frac{p_i + a_i \rho_i^2}{\gamma_i - 1}
      + \left( \frac{1 - b_i \rho_i}{\gamma - 1} - 1 \right) 2 a_i \rho_i 
$$












### Tait's equation of state

For a material obeying the Tait EOS, the pressure is defined as
$$
    p_i = (\gamma_i - 1) \rho_i e_i - \gamma_i (b_i - a_i)
$$

where $a_i$, $b_i$, and $\gamma_i$ are material properties. In M.G. form, we have
$$
    \Gamma_i = \gamma_i 
$$

$$
    \Pi_i = \gamma_i (b_i - a_i) 
$$

and

$$
    \delta_i = 0 
$$










### Stiffened gas

For a material obeying the stiffened EOS, the pressure is defined as

$$
    p_i = (\gamma_i - 1) \rho_i e_i - \gamma_i a_i 
$$

where $a_i$ and $\gamma_i$ are material properties, and

$$
    \Gamma_i = \gamma_i 
$$

$$
    \Pi_i = \gamma_i a_i 
$$

$$
    \delta_i = 0 
$$















### Jones Wilkins Lee (JWL)

The more complicated JWL EOS is often used to define energetic materials, and has a reference pressure given by

$$
    p_{ref,i} = A_i e^{-\frac{R_{1,i} \rho_{0,i}}{\rho_i}}
        + B_i e^{-\frac{R_{2,i} \rho_{0,i}}{\rho_i}}
$$

$A_i$, $B_i$, $R_{1,i}$, $R_{2,i}$, $\rho_{0,i}$, and $\Gamma_{0,i}$ are the material properties. The functions of the EOS are given by

$$
    \Gamma_i = \Gamma_{0,i} + 1 
$$


$$
    \Pi_i = 
        \Gamma_{0,i} \rho_i 
        \left(
            \frac{A_i}{R_{1,i} \rho_{0,i}} e^{-\frac{R_{1,i} \rho_{0,i}}{\rho_i}}
          + \frac{B_i}{R_{2,i} \rho_{0,i}} e^{-\frac{R_{2,i} \rho_{0,i}}{\rho_i}}
          + e_{0,i}
        \right)
      - p_{ref,i} 
$$

and

$$
    \delta_i =& \\
          &A_i e^{-\frac{R_{1,i} \rho_{0,i}}{\rho_i}}
            \left[
                \Gamma_{0,i}
                \left(
                    \frac{1}{R_{1,i} \rho_{0,i}} 
                  + \frac{1}{\rho_i}
                \right)
              - \frac{R_{1,i} \rho_{0,i}}{\rho_i^2}
            \right] \frac{1}{\Gamma_{0,i}} \\
          +&B_i e^{-\frac{R_{2,i} \rho_{0,i}}{\rho_i}}
            \left[
                \Gamma_{0,i}
                \left(
                    \frac{1}{R_{2,i} \rho_{0,i}} 
                  + \frac{1}{\rho_i}
                \right)
              - \frac{R_{2,i} \rho_{0,i}}{\rho_i^2}
            \right] \frac{1}{\Gamma_{0,i}} \\
          +&e_{0,i} 
$$

$e_0$ is anther material parameter which denotes a reference energy state.















### Cochran Chan

The Cochran Chan EOS can be used to describe solid material, and has a reference pressure given by
$$
    p_{ref,i} = A_i \left(\frac{\rho_{0,i}}{\rho_i}\right)^{1 - \mathcal{E}_{1,i}}
        - B_i \left(\frac{\rho_{0,i}}{\rho_i}\right)^{1 - \mathcal{E}_{2,i}}
$$

$A_i$, $B_i$, $\mathcal{E}_{1,i}$, $\mathcal{E}_{2,i}$, $\rho_{0,i}$, and $\Gamma_{0,i}$ are the material properties. The functions of the EOS are given by

$$
    \Gamma_i = \Gamma_{0,i} + 1
$$

$$
    \Pi_i = 
        \Gamma_{0,i} \rho_i 
        \left(
          - \frac{A_i}{(\mathcal{E}_{1,i} - 1) \rho_{0,i}} 
            \left(\frac{\rho_{0,i}}{\rho_i}\right)^{1 - \mathcal{E}_{1,i}}
          + \frac{B_i}{(\mathcal{E}_{2,i} - 1) \rho_{0,i}} 
            \left(\frac{\rho_{0,i}}{\rho_i}\right)^{1 - \mathcal{E}_{2,i}}
          + e_{0,i}
        \right)
      - p_{ref,i}
$$


$$
    \delta_i =& \\
          &\frac{A_i}{\mathcal{E}_{1,i}}
            \left[
                \mathcal{E}_{1,i} 
                \left(\frac{\rho_{0,i}}{\rho_i}\right)^{-\mathcal{E}_{1,i}}
                \frac{\mathcal{E}_{1,i} - \Gamma_{0,i} - 1}{\rho_i}
              + \frac{\Gamma_{0,i}}{\rho_{0,i}}
            \right] \frac{1}{\Gamma_{0,i}} \\
          +&\frac{B_i}{\mathcal{E}_{2,i}}
            \left[
                \mathcal{E}_{2,i} 
                \left(\frac{\rho_{0,i}}{\rho_i}\right)^{-\mathcal{E}_{2,i}}
                \frac{\mathcal{E}_{2,i} - \Gamma_{0,i} - 1}{\rho_i}
              + \frac{\Gamma_{0,i}}{\rho_{0,i}}
            \right] \frac{1}{\Gamma_{0,i}} \\
          +&e_{0,i} 
$$

where again, $e_0$ is a material parameter which denotes a reference energy state.





