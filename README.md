# ModelCloud
Routine to compute the radiation profiles, and the radiation/gravity forces, for the idealised model presented in Appendix C of Menon et al 2022b (M22). This routine solves an ordinary differential equation (ODE; Equation C5 in M22) for the solution of the radiation energy and opacity profiles. 
The routine allows for the following opacity laws: Constant opacity, Semenov opacities (Semenov 2003), and power-law approximation to the Semenov opacities. 
This repository uses some borrowed code from the ODE solution routine provided in the repository of the Quokka radiation hydrodynamics code (see Refs below), and is based on the equations described in Skinner & Ostriker 2013.
It also uses a modified version of the publicly available routine that calculates the Rosseland opacity given an input temperature and density of gas.

User parameters:
1. $\Sigma$: the (total) mass surface density of the cloud
2. $\alpha$: the power-law exponent for the radial profile of density, i.e. $\rho \propto r^{-\alpha}$
3. $\epsilon_{*}$: the fraction of mass in the stars/sources

There are also some fixed parameters, such as the total mass in the cloud, size of source etc, that are _hardcoded_. 

## Usage
If using the Semenov opacities, one would have to run `gmake` for compiling the Semenov opacity calculation routine. This can be done with - 
```
make all
```
Please check the Makefile and compiler specifications if you have issues with this. The routine to run the model can first be made an executable with 
```
chmod +x IdealCloud.py
```
and then run as
```
IdealCloud -Sigma sigmavalue -alpha alphavalue
```
where
1. `sigmavalue`: the target $\Sigma$ of the cloud in $M_{\odot} \, \mathrm{pc}^{-2}$. Default is $3.2 \times 10^3 \, M_{\odot} \, \mathrm{pc}^{-2}$
2. `alphavalue`: the power-law exponent for the density profile. Options are $\alpha = 0,1,2$. Default is $\alpha=0$, i.e. constant density.

The routine should then solve the ODE's, plot the profiles of the radiation quantities and the Eddington ratio, and save these in a text file. If the ODE has math/overflow issues, try a larger value of $\sigma_*$, parameterised as `sigma_star` in the routine.

References:
1. Semenov 2003: https://ui.adsabs.harvard.edu/abs/2003A%26A...410..611S/abstract
2. Menon et al 2022b: TBD
3. Quokka RHD code: https://github.com/BenWibking/quokka
4. Skinner & Ostriker 2013: https://ui.adsabs.harvard.edu/abs/2013ApJS..206...21S/abstract
5. Semenov opacity routine: https://www2.mpia-hd.mpg.de/~semenov/Opacities/opacities.html

