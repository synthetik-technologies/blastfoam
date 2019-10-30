<!--- python3 -m readme2tex --nocdn --pngtrick --htmlize --output README.md INPUT.md && open README.md.html -->


<p align="center">
  <img src="media/synthetik-logo.png" width="350" title="hover text">
</p>


# blastFoam

blastFoam is a solver for multi-component compressible flow with application to high-explosive detonation, explosive safety and airblast. blastFoam is developed by [Synthetik Applied Technologies](https://www.synthetik-technologies.com). 



## Disclaimer

This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM<img alt="$\textregistered$" src="svgs/6abda71802c3922eebfcf1b67d5169b2.png" align="middle" width="16.438455000000005pt" height="22.831379999999992pt"/>  and OpenCFD<img alt="$\textregistered$" src="svgs/6abda71802c3922eebfcf1b67d5169b2.png" align="middle" width="16.438455000000005pt" height="22.831379999999992pt"/> trade marks.



## What is blastFoam?

blastFoam is a solver for multi-component compressible flow with application to high-explosive detonation, explosive safety and airblast. The blastFoam solver uses OpenFOAM technology, and is in no way approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM<img alt="$\textregistered$" src="svgs/6abda71802c3922eebfcf1b67d5169b2.png" align="middle" width="16.438455000000005pt" height="22.831379999999992pt"/>  and OpenCFD<img alt="$\textregistered$" src="svgs/6abda71802c3922eebfcf1b67d5169b2.png" align="middle" width="16.438455000000005pt" height="22.831379999999992pt"/> trade marks (see Disclaimer above).



## How to use blastFoam

Several validation and tutorial cases are included in the repository, and are documented in the [blastFoam User Guide](blastFoam_User_Guide.pdf).



## Installation

1. Install OpenFOAM from [source](https://github.com/OpenFOAM/OpenFOAM-7) or via binary/package installation from [here](https://openfoam.org/version/7/) (blastFoam works with OF version 7).

2. Clone the blastFoam repository

```bash
git clone https://github.com/synthetik-technologies/blastfoam.git
```

3. Run the ./Allwmake script to compile and install blastFoam

```bash
cd blastFoam
./Allwmake
```

4. Test your installation by running the tutorial and validation cases



## Questions and Availability
If you find any bugs, please let us know in the issues section of the repository. If you want to get in touch: info@synthetik-technologies.com 

blastFoam is also available on the Texas Advanced Computing Center https://www.tacc.utexas.edu (TACC)




## Citation

If you use this code for your work or research, please use this citation:

```
J. Heylmun, P. Vonk, and T. Brewer, blastFoam: An OpenFOAM Solver for Compressible Multi-Fluid Flow with Application to High-Explosive Detonation. Synthetik Applied Technologies, LLC., 2019.
```
BiBTex:
```
@software{heylmun_blastfoam:_2019,
	title = {{blastFoam}: A Solver for Compressible Multi-Fluid Flow with Application to High-Explosive Detonation},
	url = {https://github.com/synthetik-technologies/blastfoam},
	publisher = {Synthetik Applied Technologies, {LLC}.},
	author = {Heylmun, Jeffrey and Vonk, Peter and Brewer, Timothy},
	date = {2019-10-22}
}
```

To cite the blastFoam User Guide:
```
J. Heylmun, P. Vonk, and T. Brewer, “blastFoam User Guide.” Synthetik Applied Technologies, LLC., 30-Oct-2019.
```
BiBTex:
```
@misc{heylmun_blastfoamguide_2019,
	title = {{blastFoam} {User} {Guide}},
	url = {https://github.com/synthetik-technologies/blastfoam},
	language = {English},
	publisher = {Synthetik Applied Technologies, LLC.},
	author = {Heylmun, Jeffrey and Vonk, Peter and Brewer, Timothy},
	month = oct,
	year = {2019}
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
- Van der Waals



## Verification and Validation

blastFoam has been validated against known solutions to standard gas dynamics problems, and against data from physical tests. Validation cases are included with example/tutorial cases as part of the solver source code. 





### Validation/Example: Internal Detonation 

Reference:
```
Joachim, Charles E., Gordon W. McMahon, Christo V. Lunderman, and Sharon B. Garner. 1999. “Airblast Effects Research: Small-Scale Experiments and Calculations.” DTIC Document.
```


Validation against experimental and simulated (CTH) data as reported in Joachim et. al.; see the paper for an explanation of scaling.

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

<p align="center"><img alt="$$&#10;    \partial_t \mathbf{U} + \nabla \cdot \mathbf{F} = \mathbf{S}&#10;$$" src="svgs/10b5fb69900ce3e0ba80d06da2efff49.png" align="middle" width="119.034795pt" height="13.881251999999998pt"/></p>

<p align="center"><img alt="$$&#10;    \mathbf{U} = &#10;        \left( \begin{array}{c}&#10;           \alpha_1 \\&#10;           \alpha_1 \rho_1 \\&#10;           \alpha_2 \rho2 \\&#10;           \rho \mathbf{u} \\&#10;           \rho E&#10;        \end{array} \right)&#10;    \mathbf{F} = &#10;        \left( \begin{array}{c}&#10;           \alpha_1 \mathbf{u} \\&#10;           \alpha_1 \rho_1 \mathbf{u} \\&#10;           \alpha_2 \rho_2 \mathbf{u} \\&#10;           \rho \mathbf{u} \otimes \mathbf{u} + p \mathbf{I}\\&#10;           (\rho E + p) \mathbf{u}&#10;        \end{array} \right)&#10;    \mathbf{S} = &#10;        \left( \begin{array}{c}&#10;           \alpha_1 \nabla \cdot \mathbf{u} \\&#10;           0 \\&#10;           0 \\&#10;           0 \\&#10;           0&#10;        \end{array} \right)&#10;$$" src="svgs/e17a3de7478e3ea4b1124d2aaf70c7cb.png" align="middle" width="417.48629999999997pt" height="98.63106pt"/></p>

where <img alt="$\rho$" src="svgs/6dec54c48a0438a5fcde6053bdb9d712.png" align="middle" width="8.498985000000003pt" height="14.155350000000013pt"/> is the mixture density, <img alt="$\mathbf{u}$" src="svgs/129c5b884ff47d80be4d6261a476e9f1.png" align="middle" width="10.502250000000002pt" height="14.61206999999998pt"/> the mixture velocity, <img alt="$E$" src="svgs/84df98c65d88c6adf15d4645ffa25e47.png" align="middle" width="13.082190000000004pt" height="22.46574pt"/> the total energy, <img alt="$p$" src="svgs/2ec6e630f199f589a2402fdf3e0289d5.png" align="middle" width="8.270625000000004pt" height="14.155350000000013pt"/> the pressure, and <img alt="$rho_i$" src="svgs/b921dbff3fd97c42704607ee82ee7cc8.png" align="middle" width="29.96301pt" height="22.831379999999992pt"/> and <img alt="$\alpha_i$" src="svgs/67e7dd600dde6ca2d15b4df76a96110b.png" align="middle" width="15.166635000000003pt" height="14.155350000000013pt"/> are the density and volume fraction of each phase. 





<p align="center"><img alt="$$&#10;    \alpha_2 = 1 - \alpha_1&#10;$$" src="svgs/7095db6a8ced11129e3a9e734856f850.png" align="middle" width="85.186365pt" height="13.059337499999998pt"/></p>


<p align="center"><img alt="$$&#10;    \rho = \sum_i \alpha_i \rho_i&#10;$$" src="svgs/f0a7c0cae2506902c7bdae655786b0c4.png" align="middle" width="86.038755pt" height="36.655409999999996pt"/></p>


<p align="center"><img alt="$$&#10;    \rho E = \rho e + \frac{1}{2}\rho |\mathbf{u}|^2&#10;$$" src="svgs/17c1f76bf866c00914f7f2c4224ef392.png" align="middle" width="126.59344499999999pt" height="32.9901pt"/></p>

<p align="center"><img alt="$$&#10;    \rho e = \sum_i \alpha_i \rho_i e_i&#10;$$" src="svgs/f3a115a9e214cc2aaa59556c40aa818f.png" align="middle" width="106.819845pt" height="36.655409999999996pt"/></p>

The pressure will be defined using a specified equation of state where the mixture internal energy, densities, and volume fraction are used to calculate the total pressure. The equations of state will be used in the Mie-Gruneisen form.

<p align="center"><img alt="$$&#10;    p_i(\rho_i, e_i, \rho) = (\Gamma(\rho_i) - 1) \rho_i e_i - \Pi(\rho_i)&#10;$$" src="svgs/5d9c11f4f060e305ea1696832228d0d0.png" align="middle" width="263.05785pt" height="16.438356pt"/></p>







## Equations of State 


The mixture pressure is defined using the Mie-Gruneisen from using
<p align="center"><img alt="$$&#10;    p = \frac{\rho e}{\sum_i \alpha_i \xi_i} - \frac{\sum_i \alpha_i \xi_i \Pi_i}{\sum_i \alpha_i \xi_i}&#10;$$" src="svgs/a4b9c8c1b20f00ca4b636b1c8d0a792a.png" align="middle" width="182.4339pt" height="39.878685pt"/></p> 

where 

<p align="center"><img alt="$$&#10;    \xi_i(\rho_i) = \frac{1}{\Gamma_i - 1}&#10;$$" src="svgs/74e931d77f7803a5d2fbf6e2c4482778.png" align="middle" width="107.36912999999998pt" height="35.455859999999994pt"/></p>

and <img alt="$\Pi_i$" src="svgs/039d0f21249b35a3d8cf7be5cdb855fc.png" align="middle" width="16.979820000000004pt" height="22.46574pt"/> is dependent on the equation of state.

The speed of sound within a give phase is give by

<p align="center"><img alt="$$&#10;    c_i = \sqrt{\frac{\sum_i y_i \xi_i c_i^2}{\sum_i \xi_i}}&#10;$$" src="svgs/e8f5da9a7fdd1ddfb8264a0627e4c209.png" align="middle" width="121.13705999999999pt" height="49.31553pt"/></p>

with

<p align="center"><img alt="$$&#10;    y_i = \frac{\alpha_i \rho_i}{\rho}&#10;$$" src="svgs/10d0eca7b36d3bf77906b7aad1b3977a.png" align="middle" width="67.38253499999999pt" height="32.670494999999995pt"/></p>



<p align="center"><img alt="$$&#10;    c_i^2 = \frac{h_i - \delta_i}{\xi_i}&#10;$$" src="svgs/2073832399f0bd906131d86a3c7d06d1.png" align="middle" width="86.19237pt" height="37.00851pt"/></p>



<p align="center"><img alt="$$&#10;    h_i = \frac{\Gamma_i p + \Pi_i}{(\Gamma_i - 1)\rho_i}&#10;$$" src="svgs/df2308ed019830f3a532ef0ca0c2c5cc.png" align="middle" width="109.64843999999998pt" height="37.738965pt"/></p>

and <img alt="$\delta_i$" src="svgs/8ed540494c9132490e49ade18ca8d273.png" align="middle" width="11.956890000000003pt" height="22.831379999999992pt"/> is again dependent on the equation of state.












### Generalized van der Waals gas

For a gas described by the generalized van der Waals equation of state, the pressure is defined as
<p align="center"><img alt="$$&#10;    p_i = \frac{\gamma_i - 1}{1 - b_i \rho_i}(\rho_i e_i + a_i \rho_i^2) - (a_i \rho_i^2 + c_i) &#10;$$" src="svgs/12e2ad5186082905d82b059b18242141.png" align="middle" width="280.01655pt" height="36.186479999999996pt"/></p>

where <img alt="$a_i$" src="svgs/65ed4b231dcf18a70bae40e50d48c9c0.png" align="middle" width="13.340085000000004pt" height="14.155350000000013pt"/>, <img alt="$b_i$" src="svgs/d3aa71141bc89a24937c86ec1d350a7c.png" align="middle" width="11.705760000000003pt" height="22.831379999999992pt"/>, <img alt="$c_i$" src="svgs/3bc6fc8b86b6c61889f4e572c7546b8e.png" align="middle" width="11.764830000000003pt" height="14.155350000000013pt"/>, and <img alt="$\gamma_i$" src="svgs/9925e4d3a0486c8d876c6c8eec9e256d.png" align="middle" width="13.161555000000003pt" height="14.155350000000013pt"/> are material parameters. Writing the above equation in M.G. form, we obtain

<p align="center"><img alt="$$&#10;    \Gamma_i = \frac{\gamma_i - 1}{1 - b_i \rho_i} + 1 &#10;$$" src="svgs/78b3bdddaf7ba69d93d007ee3e7f378d.png" align="middle" width="124.72960499999999pt" height="36.186479999999996pt"/></p>

<p align="center"><img alt="$$&#10;    \Pi_i = &#10;        \left[1 - \frac{\gamma_i - 1}{1 - b_i \rho_i}\right] a_i \rho_i^2 &#10;      + \left[\frac{\gamma_i - 1}{1 - b_i \rho_i} + 1\right] c_i &#10;$$" src="svgs/7ab1654988023cd894bf926eda793c96.png" align="middle" width="315.92384999999996pt" height="39.45249pt"/></p>

and

<p align="center"><img alt="$$&#10;    \delta_i = &#10;      - b_i \frac{p_i + a_i \rho_i^2}{\gamma_i - 1}&#10;      + \left( \frac{1 - b_i \rho_i}{\gamma - 1} - 1 \right) 2 a_i \rho_i &#10;$$" src="svgs/9f0af7e53e3ca2854fd8d65af0d90b5d.png" align="middle" width="297.4521pt" height="40.118265pt"/></p>












### Tait's equation of state

For a material obeying the Tait EOS, the pressure is defined as
<p align="center"><img alt="$$&#10;    p_i = (\gamma_i - 1) \rho_i e_i - \gamma_i (b_i - a_i)&#10;$$" src="svgs/82847845b0e03f582d590ea228f5d359.png" align="middle" width="211.48049999999998pt" height="16.438356pt"/></p>

where <img alt="$a_i$" src="svgs/65ed4b231dcf18a70bae40e50d48c9c0.png" align="middle" width="13.340085000000004pt" height="14.155350000000013pt"/>, <img alt="$b_i$" src="svgs/d3aa71141bc89a24937c86ec1d350a7c.png" align="middle" width="11.705760000000003pt" height="22.831379999999992pt"/>, and <img alt="$\gamma_i$" src="svgs/9925e4d3a0486c8d876c6c8eec9e256d.png" align="middle" width="13.161555000000003pt" height="14.155350000000013pt"/> are material properties. In M.G. form, we have
<p align="center"><img alt="$$&#10;    \Gamma_i = \gamma_i &#10;$$" src="svgs/5a10fc7dadbf5346b033fd11058a609e.png" align="middle" width="50.825939999999996pt" height="14.429217pt"/></p>

<p align="center"><img alt="$$&#10;    \Pi_i = \gamma_i (b_i - a_i) &#10;$$" src="svgs/13d62767165958a1e91b642f9dd36ec3.png" align="middle" width="113.26887pt" height="16.438356pt"/></p>

and

<p align="center"><img alt="$$&#10;    \delta_i = 0 &#10;$$" src="svgs/62e8a915009a0f0acc7a12cadb16551b.png" align="middle" width="42.91551pt" height="13.881251999999998pt"/></p>










### Stiffened gas

For a material obeying the stiffened EOS, the pressure is defined as

<p align="center"><img alt="$$&#10;    p_i = (\gamma_i - 1) \rho_i e_i - \gamma_i a_i &#10;$$" src="svgs/91090e3a49d9547c75589fd44068bd6b.png" align="middle" width="165.2541pt" height="16.438356pt"/></p>

where <img alt="$a_i$" src="svgs/65ed4b231dcf18a70bae40e50d48c9c0.png" align="middle" width="13.340085000000004pt" height="14.155350000000013pt"/> and <img alt="$\gamma_i$" src="svgs/9925e4d3a0486c8d876c6c8eec9e256d.png" align="middle" width="13.161555000000003pt" height="14.155350000000013pt"/> are material properties, and

<p align="center"><img alt="$$&#10;    \Gamma_i = \gamma_i &#10;$$" src="svgs/5a10fc7dadbf5346b033fd11058a609e.png" align="middle" width="50.825939999999996pt" height="14.429217pt"/></p>

<p align="center"><img alt="$$&#10;    \Pi_i = \gamma_i a_i &#10;$$" src="svgs/3cc94a60a6f036e7ab481d51bdcf2495.png" align="middle" width="67.0428pt" height="14.429217pt"/></p>

<p align="center"><img alt="$$&#10;    \delta_i = 0 &#10;$$" src="svgs/62e8a915009a0f0acc7a12cadb16551b.png" align="middle" width="42.91551pt" height="13.881251999999998pt"/></p>















### Jones Wilkins Lee (JWL)

The more complicated JWL EOS is often used to define energetic materials, and has a reference pressure given by

<p align="center"><img alt="$$&#10;    p_{ref,i} = A_i e^{-\frac{R_{1,i} \rho_{0,i}}{\rho_i}}&#10;        + B_i e^{-\frac{R_{2,i} \rho_{0,i}}{\rho_i}}&#10;$$" src="svgs/93a9ffd37dc4fcbea984096ae3e6f704.png" align="middle" width="244.8204pt" height="26.656575pt"/></p>

<img alt="$A_i$" src="svgs/4ebf880807deff5796460f39aea46f80.png" align="middle" width="16.979820000000004pt" height="22.46574pt"/>, <img alt="$B_i$" src="svgs/72f4aab7f49593ada1f6b406b90a8a94.png" align="middle" width="17.119575000000005pt" height="22.46574pt"/>, <img alt="$R_{1,i}$" src="svgs/63ca638312fae416c900385b2fa19d49.png" align="middle" width="27.589155000000005pt" height="22.46574pt"/>, <img alt="$R_{2,i}$" src="svgs/b5d4de368d7931bfe1f4fa2aaf82f2ac.png" align="middle" width="27.589155000000005pt" height="22.46574pt"/>, <img alt="$\rho_{0,i}$" src="svgs/32f611e8510475eae12d6fea66e3d5ea.png" align="middle" width="23.606550000000002pt" height="14.155350000000013pt"/>, and <img alt="$\Gamma_{0,i}$" src="svgs/6c14aae54ace497cf43fcfc933edfcc5.png" align="middle" width="25.381620000000005pt" height="22.46574pt"/> are the material properties. The functions of the EOS are given by

<p align="center"><img alt="$$&#10;    \Gamma_i = \Gamma_{0,i} + 1 &#10;$$" src="svgs/4831c30c7cf835f3433ca81f63bcec8b.png" align="middle" width="92.17824pt" height="15.93603pt"/></p>


<p align="center"><img alt="$$&#10;    \Pi_i = &#10;        \Gamma_{0,i} \rho_i &#10;        \left(&#10;            \frac{A_i}{R_{1,i} \rho_{0,i}} e^{-\frac{R_{1,i} \rho_{0,i}}{\rho_i}}&#10;          + \frac{B_i}{R_{2,i} \rho_{0,i}} e^{-\frac{R_{2,i} \rho_{0,i}}{\rho_i}}&#10;          + e_{0,i}&#10;        \right)&#10;      - p_{ref,i} &#10;$$" src="svgs/f05bd689f4e7e0427425e8ede6a7c546.png" align="middle" width="473.3025pt" height="39.81482999999999pt"/></p>

and

<p align="center"><img alt="$$&#10;    \delta_i =&amp; \\&#10;          &amp;A_i e^{-\frac{R_{1,i} \rho_{0,i}}{\rho_i}}&#10;            \left[&#10;                \Gamma_{0,i}&#10;                \left(&#10;                    \frac{1}{R_{1,i} \rho_{0,i}} &#10;                  + \frac{1}{\rho_i}&#10;                \right)&#10;              - \frac{R_{1,i} \rho_{0,i}}{\rho_i^2}&#10;            \right] \frac{1}{\Gamma_{0,i}} \\&#10;          +&amp;B_i e^{-\frac{R_{2,i} \rho_{0,i}}{\rho_i}}&#10;            \left[&#10;                \Gamma_{0,i}&#10;                \left(&#10;                    \frac{1}{R_{2,i} \rho_{0,i}} &#10;                  + \frac{1}{\rho_i}&#10;                \right)&#10;              - \frac{R_{2,i} \rho_{0,i}}{\rho_i^2}&#10;            \right] \frac{1}{\Gamma_{0,i}} \\&#10;          +&amp;e_{0,i} &#10;$$" src="svgs/abd072e5bd492908795f171d2f8d7a17.png" align="middle" width="806.88465pt" height="39.81482999999999pt"/></p>

<img alt="$e_0$" src="svgs/824a0be2acf2b955fb812cf516864cf4.png" align="middle" width="14.206830000000004pt" height="14.155350000000013pt"/> is anther material parameter which denotes a reference energy state.















### Cochran Chan

The Cochran Chan EOS can be used to describe solid material, and has a reference pressure given by
<p align="center"><img alt="$$&#10;    p_{ref,i} = A_i \left(\frac{\rho_{0,i}}{\rho_i}\right)^{1 - \mathcal{E}_{1,i}}&#10;        - B_i \left(\frac{\rho_{0,i}}{\rho_i}\right)^{1 - \mathcal{E}_{2,i}}&#10;$$" src="svgs/2ca462d36387f32e82ad4779afbb4b2e.png" align="middle" width="302.6826pt" height="43.251615pt"/></p>

<img alt="$A_i$" src="svgs/4ebf880807deff5796460f39aea46f80.png" align="middle" width="16.979820000000004pt" height="22.46574pt"/>, <img alt="$B_i$" src="svgs/72f4aab7f49593ada1f6b406b90a8a94.png" align="middle" width="17.119575000000005pt" height="22.46574pt"/>, <img alt="$\mathcal{E}_{1,i}$" src="svgs/2c2751f493c433dbaf206ac3d0e699a5.png" align="middle" width="23.78343pt" height="22.46574pt"/>, <img alt="$\mathcal{E}_{2,i}$" src="svgs/12e2967cc57a8c6605a2facfb3c45551.png" align="middle" width="23.78343pt" height="22.46574pt"/>, <img alt="$\rho_{0,i}$" src="svgs/32f611e8510475eae12d6fea66e3d5ea.png" align="middle" width="23.606550000000002pt" height="14.155350000000013pt"/>, and <img alt="$\Gamma_{0,i}$" src="svgs/6c14aae54ace497cf43fcfc933edfcc5.png" align="middle" width="25.381620000000005pt" height="22.46574pt"/> are the material properties. The functions of the EOS are given by

<p align="center"><img alt="$$&#10;    \Gamma_i = \Gamma_{0,i} + 1&#10;$$" src="svgs/fcaf18d66765c09f623992e25221d14a.png" align="middle" width="92.17824pt" height="15.93603pt"/></p>

<p align="center"><img alt="$$&#10;    \Pi_i = &#10;        \Gamma_{0,i} \rho_i &#10;        \left(&#10;          - \frac{A_i}{(\mathcal{E}_{1,i} - 1) \rho_{0,i}} &#10;            \left(\frac{\rho_{0,i}}{\rho_i}\right)^{1 - \mathcal{E}_{1,i}}&#10;          + \frac{B_i}{(\mathcal{E}_{2,i} - 1) \rho_{0,i}} &#10;            \left(\frac{\rho_{0,i}}{\rho_i}\right)^{1 - \mathcal{E}_{2,i}}&#10;          + e_{0,i}&#10;        \right)&#10;      - p_{ref,i}&#10;$$" src="svgs/f32856152e6140939fa5ed3d1fffb8fb.png" align="middle" width="619.2053999999999pt" height="49.31553pt"/></p>


<p align="center"><img alt="$$&#10;    \delta_i =&amp; \\&#10;          &amp;\frac{A_i}{\mathcal{E}_{1,i}}&#10;            \left[&#10;                \mathcal{E}_{1,i} &#10;                \left(\frac{\rho_{0,i}}{\rho_i}\right)^{-\mathcal{E}_{1,i}}&#10;                \frac{\mathcal{E}_{1,i} - \Gamma_{0,i} - 1}{\rho_i}&#10;              + \frac{\Gamma_{0,i}}{\rho_{0,i}}&#10;            \right] \frac{1}{\Gamma_{0,i}} \\&#10;          +&amp;\frac{B_i}{\mathcal{E}_{2,i}}&#10;            \left[&#10;                \mathcal{E}_{2,i} &#10;                \left(\frac{\rho_{0,i}}{\rho_i}\right)^{-\mathcal{E}_{2,i}}&#10;                \frac{\mathcal{E}_{2,i} - \Gamma_{0,i} - 1}{\rho_i}&#10;              + \frac{\Gamma_{0,i}}{\rho_{0,i}}&#10;            \right] \frac{1}{\Gamma_{0,i}} \\&#10;          +&amp;e_{0,i} &#10;$$" src="svgs/a2bfd09944f0d0622e5746946fdf45a0.png" align="middle" width="787.3585499999999pt" height="49.31553pt"/></p>

where again, <img alt="$e_0$" src="svgs/824a0be2acf2b955fb812cf516864cf4.png" align="middle" width="14.206830000000004pt" height="14.155350000000013pt"/> is a material parameter which denotes a reference energy state.





