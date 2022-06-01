#!/usr/bin/env python
# -*- coding: utf-8 -*-
# written by Shyam Menon, 2022
"""
Routine to compute the steady state profile of radiation energy density for a spherically symmetric cloud with a luminous source at the centre.
Based on the model presented in Skinner & Ostriker (2013) in their thin-shell expansion test.

Written by Ben Wibking (2021; https://github.com/BenWibking/quokka/blob/development/extern/dust_shell/solution.py)
Extended to a temperature-dependent Semenov opacity by Shyam Menon (2022)

See Appendix
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import root_scalar
from math import pi, sqrt, exp, log
import mpmath as mp

from sympy import false
import argparse

#Parse user arguments
parser = argparse.ArgumentParser(description='Options for idealised spherical cloud model.')
parser.add_argument("-Sigma","--Sigma",type=float,default=3.2e3,help="The target mass surface density of the cloud in Solar Mass/parsec^2. Default is 3.2 x 10^3.")
parser.add_argument("-alpha","--alpha",type=float,default=0,help="The power-law index such that density ~ r^{-alpha}. Default is zero, i.e. constant density sphere.")
parser.add_argument("-epsilon","--epsilon",type=float,default=0.7,help="The fraction of mass in stars. Default is 70%.")
parser.add_argument("-show","--show",action='store_true',help="Flag to show the plots interactively before saving")
args = parser.parse_args()

#Set user parameters
SigmaSolPc = args.Sigma # Surface density in Msol/pc^2
rho_alpha = args.alpha #Density profile power law slope
epsilon = args.epsilon # SFE
show = args.show
print("Input Values are: Sigma = {} Msol/pc^2 \t alpha = {} \t epsilon = {}".format(SigmaSolPc,rho_alpha,epsilon))

################################################################################
#First define some physical constants
a_rad = 7.5646e-15      # erg cm^-3 K^-4
c = 2.99792458e10       # cm s^-1    
Msun = 2.0e33           # g
parsec_in_cm = 3.086e18  # cm
G_GravityConstant = 6.67428e-08

#Some fixed problem parameters
specific_luminosity = 1.7e3             # erg s^-1 g^-1
GMC_mass = 1.0e6 * Msun                 # g
M_g = (1 - epsilon) * GMC_mass      # Gas mass in cloud in g
M_star = epsilon * GMC_mass      # Stellar mass in cloud in g
L_star = M_star * specific_luminosity  # erg s^-1
Sigma = SigmaSolPc * Msun/(parsec_in_cm**2) #Mass surface density of cloud

#Derived radius to match target Sigma
R_g = np.sqrt(GMC_mass/(np.pi*Sigma)) # Radius of cloud in cm

#radial resolution and grid
Nr = 2048
sigma_star = R_g/(Nr) * 64.0
rmin = sigma_star
rmax = R_g
# Density assuming uniform sphere
rho_0 = M_g / ((4. / 3.) * np.pi * R_g * R_g * R_g)  # g cm^-3

#Dust opacity law: options are fixed, semenov and plt - power law 
kappalaw = "semenov"
#Value if using Uniform opacity
kappa0 = 5.0           # specific opacity [cm^2 g^-1]

#Construct the model name
def find_exp(number) -> int:
    base10 = np.log10(abs(number))
    return abs(np.floor(base10))
Sig_string = str(int(find_exp(SigmaSolPc)))
Dens_string = str(int(rho_alpha))
Model_Name = "Sigma{}_Dens{}_eps{}".format(Sig_string,Dens_string,epsilon)
print("Model Name = {}".format(Model_Name))
print("")




################################################################################

def call_opacity(temp,rho):
    import subprocess
    cmd = "./rosseland {:2.3f} {:2.3e}".format(float(temp),float(rho))
    proc = subprocess.Popen([cmd],stdout=subprocess.PIPE,shell=True)
    (out, err) = proc.communicate()
    try:
        opac = float(out.split(b'output kappa ')[1].split(b'=')[1].split(b'\n')[0])
    except:
        print(cmd)
    return opac


def get_kappa(r,f):
    #r is in dimensional units
    if(kappalaw.lower() == "uniform"):
        return kappa0
    #Power-law temperature dependence
    elif(kappalaw.lower() == "plt"):
        Er = Frad(r)/(c * f)
        Tr = (Er/a_rad)**0.25
        return 10**(-1.5) * (Tr/10.0)**2
    #Semenov opacities
    elif(kappalaw.lower() == "semenov"):
        Er = Frad(r)/(c * f)
        Tr = (Er/a_rad)**0.25
        Tr = min(max(5.0,Tr),1.e4)
        dens = rho(r)
        dens = min(max(dens,2.e-18),2.e-7)
        return call_opacity(Tr,dens)

# solution for radiation flux F in cgs
def Frad(r):
    from scipy import special
    # r is in cgs!
    normfac = L_star/(4.0*np.pi*r*r)
    term1 = special.erf(r/(np.sqrt(2.0)*sigma_star))
    term2 = 2.0*r/(np.sqrt(2.0*np.pi)*sigma_star) * \
        np.exp(-0.5 * (r/sigma_star)**2)
    return normfac*(term1 - term2)

def dlnF_dr(r):
    from scipy import special
    # r is in cgs!
    # compute d ln Frad / dr
    return (-2/r) + (2*(r/sigma_star)**2) / \
        (-2*r + np.exp(0.5*(r/sigma_star)**2) * (np.sqrt(2*pi)
         * sigma_star) * special.erf(r/(np.sqrt(2)*sigma_star)))

def df_dr_lim(f, r_in):
    from scipy import special
    # r_in is in units of R_g
    r = r_in * R_g  # cm
    rhofunc = np.vectorize(rho)
    kappa = get_kappa(r,f)
    tau = kappa * rhofunc(r)
    G = dlnF_dr(r)
    sigma = sigma_star
    deriv = (2*r*(2*(10 - 4*sqrt(4 - 3*f**2) + 3*f**2*(-5 + 3*sqrt(4 - 3*f**2)))*r**2 +
                  2*(f**2*(6 - 9*sqrt(4 - 3*f**2)) + 4*(-1 + sqrt(4 - 3*f**2)))*sigma**2 +
                  3*f*(-8 + 9*f**2)*r*sigma**2*tau) +
             exp(r**2/(2.*sigma**2))*sqrt(2*pi)*sigma**3 *
             (8 - 8*sqrt(4 - 3*f**2) -
              3*f*(-8*r*tau + f*(4 - 6*sqrt(4 - 3*f**2) + 9*f*r*tau))) *
             special.erf(r/(sqrt(2)*sigma))) /\
        (15.*f*r*sigma**2*(2*r - exp(r**2/(2.*sigma**2))*sqrt(2*pi)*sigma *
                           special.erf(r/(sqrt(2)*sigma))))
    return -R_g*deriv

def df_dr(f, r_in):
    # r_in is dimensionless
    sqrtfac = np.sqrt(4.0 - 3.0*f*f)
    normfac = 3.0*f*sqrtfac / (5.0*sqrtfac - 8.0)
    kappa = get_kappa(r_in*R_g,f)
    r = r_in * R_g
    term1 = (dlnF_dr(r) / 3.0) * (5.0 - 2.0*sqrtfac)
    term2 = (2.0/r) * (2.0 - sqrtfac)
    term3 = rho(r) * kappa * f
    return R_g * normfac*(term1 + term2 + term3)

# compute ODE solution for f == F/cE (dimensionless)
def func_M1(t, y):
    # return dy / dt == func(t,y)
    # in this case, y == f (reduced flux)
    #               t == r (radius)
    flux = y[0]
    r = t
    if(flux <= 0.):
        print(f"r = {r}; f = {flux}")
        assert(flux > 0.)
    deriv = df_dr(flux, r)
    return [deriv]

# density profile in cgs. Uses mass conservation to get density profile
def rho(r):
    # r is in cgs!
    # compute density profile
    alpha = rho_alpha
    #Constant density
    if(alpha == 0):
        if(r<R_g):
            return rho_0
        else:
            return 1.e-10*rho_0
    if(alpha == 1):
        if(r>R_g):
            return 1.e-10*rho_0
        A = M_g/(2.*np.pi*R_g**2)
        return A/r
    if(alpha == 2):
        if(r>R_g):
            return 1.e-10*rho_0
        A = M_g/(4*np.pi*R_g)
        return A/r**2
    else:
        raise ValueError("The case for alpha = {} not considered".format(alpha))

# solve for critical point
def g(r):
    rhofunc = np.vectorize(rho)
    kappafunc = np.vectorize(get_kappa,excluded='f')
    kappa = kappafunc(r*R_g,f=2*sqrt(3)/5.)
    return 3.0*R_g*dlnF_dr(r*R_g) + (4.0/r) + 2.0*sqrt(3.0)*R_g*rhofunc(r*R_g)*kappa

def solve_branch(f0, r0, r1):
        # this equation is stiff
        sol = solve_ivp(func_M1, [r0, r1], [f0], method='BDF')
        #print(f"{sol}")
        return sol

def cumul_Mass(r,rArray):
    indices = np.where(rArray<=r)
    rhofunc = np.vectorize(rho)
    dens = rhofunc(rArray)
    y = 4*np.pi*rArray[indices]**2 * dens[indices]
    return np.trapz(y,rArray[indices])

if __name__ == '__main__':

    # print(rho_0,R_g/parsec_in_cm,M_star/Msun,Sigma)
    # r = np.linspace(0.1,20.0,100)
    # plt.plot(r,g(r))
    # plt.show(block=True)

    #Supress runtime overflow warnings
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        r_solution = np.linspace(rmin,rmax,Nr)

        try:
            root_sol = root_scalar(g, bracket=[rmin/R_g, rmax*3/R_g], method='bisect')
            print(f"critical point r_crit = {root_sol}\n")
        except ValueError:
            r = np.linspace(rmin/R_g,5*rmax/R_g,100)
            rhofunc = np.vectorize(rho)
            plt.plot(r,rhofunc(r))
            plt.axhline(0.0,ls='--',c='k',alpha=0.3)
            plt.xlabel(r'$r/R_{\mathrm{cloud}}$')
            plt.ylabel(r'$\rho$')
            plt.show(block=True)

            
            plt.plot(r,g(r))
            plt.axhline(0.0,ls='--',c='k',alpha=0.3)
            plt.ylim(-10.0,10.0)
            plt.xlabel('radius r (dimensionless)')
            plt.ylabel('function')
            plt.show(block=True)

            

        #Critical point is now a boundary condition
        r_crit = root_sol.root  # critical point at f(r_crit) = f_crit
        f_crit = 2.0*sqrt(3.0) / 5.0
        df_dr_crit = df_dr_lim(f_crit, r_crit)
        eps = 1.0e-6

        
        solve_left = True
        if solve_left:
            # solve left branch
            r_start = r_crit - eps
            r_end = rmin/R_g  # in units of R_g
            f_start = f_crit - eps * df_dr_crit
            print(f"r_crit = {r_crit}; r_start = {r_start}")
            print(f"f_crit = {f_crit}; f_start = {f_start}")
            print(f"df_dr(f_crit) = {df_dr_crit}\n")
            left_sol = solve_branch(f_start, r_start, r_end)

        #If critical radius within the max radius, then integrate to the right
        if(rmax>r_crit):
            solve_right=True

        if(solve_right):
            # solve right branch
            r_start = r_crit + eps
            r_end = 3.5  # in units of R_g
            f_start = f_crit + eps * df_dr_crit
            print(f"r_crit = {r_crit}; r_start = {r_start}")
            print(f"f_crit = {f_crit}; f_start = {f_start}")
            print(f"df_dr(f_crit) = {df_dr_crit}\n")
            right_sol = solve_branch(f_start, r_start, r_end)
            print("")




    #Plot solution and profiles of radiation quantities
    fig,axs = plt.subplots(nrows=3,figsize=(5,8),sharex='all',constrained_layout=True)

    r_over_r0 = np.concatenate((left_sol.t[::-1], right_sol.t))
    reduced_flux = np.concatenate((left_sol.y[0][::-1], right_sol.y[0]))


    #Eddington factor
    axs[0].axvline(rmin/R_g,ls='--',c='k',alpha=0.3)
    axs[0].axvline(rmax/R_g,ls='--',c='k',alpha=0.3)
    if solve_left:
        axs[0].plot(left_sol.t, left_sol.y[0], label='left branch',lw=3.0)
    if solve_right:
        axs[0].plot(right_sol.t, right_sol.y[0], label='right branch',lw=3.0)
    axs[0].axvline(rmin/R_g,ls='--',c='k',alpha=0.3)
    axs[0].axvline(rmax/R_g,ls='--',c='k',alpha=0.3)
    axs[0].scatter(r_crit, f_crit, color='black', s=20.0,zorder=2.5)
    axs[0].set_ylabel(r'$f = F_r/(cE_r)$')

    #Trad
    Frad_vec = np.vectorize(Frad)
    Frad_cgs = Frad_vec(r_over_r0 * R_g)
    Erad_cgs = Frad_cgs / (c*reduced_flux)
    Trad = (Erad_cgs/a_rad)**0.25
    axs[1].plot(r_over_r0,Trad,lw=3.0)
    axs[1].set_ylabel(r'$T_r \, (\mathrm{K})$')

    #Kappa
    kappafunc = np.vectorize(get_kappa)
    kappa = kappafunc(r_over_r0*R_g,reduced_flux)
    axs[2].plot(r_over_r0,kappa,lw=3.0)
    axs[2].set_ylabel(r'$\kappa_{\mathrm{R}} \, (\mathrm{cm}^2 \, \mathrm{g}^{-1})$')    


    axs[2].set_xlabel(r'$r/R_{\mathrm{cloud}}$')
    if(show):
        plt.show(block=True)
    fig.savefig("{}_profiles.pdf".format(Model_Name),bbox_inches='tight')
    plt.close(fig)

    #Compute Eddington Ratio

    #Compute gravitational forces
    g_sink = M_star*G_GravityConstant/(r_over_r0*R_g)**2
    g_self = np.zeros(np.size(r_over_r0))
    for i,r in enumerate(r_over_r0):
        g_self[i] = cumul_Mass(r*R_g,r_over_r0)
    g_self *= G_GravityConstant/(r_over_r0*R_g)**2
    g_total = g_self + g_sink

    #Radiation forces
    f_rad = Frad(r_over_r0*R_g)*kappa/c

    #Eddington Ratio
    fEdd = f_rad/g_total

    #Plot this
    fig,axs = plt.subplots(ncols=1)
    axs.plot(r_over_r0,fEdd,lw=3.0)
    #axs.plot(r_over_r0,g_total,lw=3.0,ls='--')
    axs.set_xlabel(r'$r/R_{\mathrm{cloud}}$')
    axs.set_ylabel(r'$f_{\mathrm{Edd}}$')
    if(show):
        plt.show(block=True)
    fig.savefig("{}_fEdd.pdf".format(Model_Name),bbox_inches='tight')
    plt.close(fig)

    #Save everything in a text file
    data = (r_over_r0,r_over_r0*R_g,reduced_flux,Frad_cgs,Trad,kappa,g_sink,g_self,g_total,f_rad,f_rad/g_total)
    header = '[0]:r/r0 \t [1]:r \t [2]:fred \t [3]:Fr \t [4]: Tr \t [5]: kappa \t [6]:g_sink \t [7]:g_self \t [8]:g_total \t [9]:f_rad \t [10]: fEdd'
    np.savetxt("{}.dat".format(Model_Name),data,delimiter='\t',header=header)


    







    


