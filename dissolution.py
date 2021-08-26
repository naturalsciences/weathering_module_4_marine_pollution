# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 14:00:51 2021

TODO: put friction velocity out of fernandez?

@author: Ludovic Lepers
"""
import math
import numpy as np
import evaporation as ev

def solubility_time(init_s, time, a = -0.1/3600):
    """
    Return the solubility as a time fonction [kg/m³]
    source : (Shen et al., 1993)

    Parameters
    ----------
    init_s : Initial solubility [kg/m³]
    time : Time from the start of the oil spill [s]
    a : Parameter of the exponential, the default is -0.1/3600.


    """
    return init_s * math.exp(a * time)



def solubility_fract_evap(init_s, frac, a = -12):
    """
    Return the solubility as a evaporated fraction fonction [kg/m³]
    source : (Shen et al., 1993)

    Parameters
    ----------
    init_s : Initial solubility [kg/m³]
    frac : Evaporated fraction []
    a : Parameter of the exponential, the default is -0.1/3600.


    """
    return init_s * math.exp(a * frac)



def diffusion_coefficient(mol_vol, wat_viscosity = 10**-3):
    """
    Return the diffusion coefficent [m²/s]
    source : (HNS-MS)

    Parameters
    ----------
    mol_vol : Molar volume of the component [mol/m³]
    wat_viscosity : Dynamic viscosity of water [Pa s]

    """
    return (13.26*10**-5)/((wat_viscosity*1000)**1.14 * (mol_vol*100**3)**0.589)/(100*100)


def molar_volume_Le_Bas(MW, organic = True):
    """
    Return the molar volume [m³/mol] from the molecular weight
    source : (HNS-MS)

    Parameters
    ----------
    MW : Molecular weight [kg/mol]
    organic : Chemical type, True if organic, the default is True.


    """
    if organic:
        return 4.9807  * (MW/1000)**0.6963
    else:
        return 2.8047 * (MW/1000)**0.651


def sherwood_slick(schmdt_nmbr, wind_speed, L, wat_kin_visc = 10**-6):
    """
    Return the average sherwood number for a slick []
    source : (HNS-MS)

    Parameters
    ----------
    schmdt_nmbr : Schmidt number []
    wind_speed : Wind speed (10 meters?) [m/s]
    L : Diameter of the slick [m]
    wat_kin_visc : Kinematic viscosity of water, the default is 10**-6 [m²/s]


    """
    re = wind_speed * L/ wat_kin_visc
    return 0.578 * schmdt_nmbr**(1/3) * math.sqrt(re)


def sherwood_droplet(schmdt_nmbr, r_velocity, L, wat_kin_visc = 10**-6):
    """
    Return the average sherwood number for a dropplet []
    source : (HNS-MS)

    Parameters
    ----------
    schmdt_nmbr : Schmidt number []
    wind_speed : Wind speed (10 meters?) [m/s]
    L : Diameter of the droplet [m]
    wat_kin_visc : Kinematic viscosity of water, the default is 10**-6 [m²/s]


    """
    re = r_velocity * L/ wat_kin_visc
    return 2+ 0.347 * schmdt_nmbr**(0.31) * re**0.62



def mass_transfer_coefficient_HNS(Sh, Dc, diameter):
    """
    return the mass transfer coefficent [m/s]
    source : (HNS-MS)

    Parameters
    ----------
    Sh : Average sherwood number []
    Dc : Diffusion coefficient at 25 °C[m²/s]
    diameter : Diameter of the slick or the droplet [m]

    """
    return (Sh*Dc) / diameter



def mass_transfer_coefficient_HNS_Fern(wind_speed, schmdt_nmbr, sct = 0.8, vnkarm = 0.47,
                                       rho_wat = 1000, rho_air = 1.2, wat_visc= 1e-3):
    """
    Compute the mass tranfer coefficient [m/s] from HNS with Fernandez.
    source : HNS MS

    Parameters
    ----------
    wind_speed : Wind speed [m/s]
    schmdt_nmbr : Schmidt number []
    sct : Turbulent schmidt number is 0.8
    vnkarm : Von karman constant = 0.47
    rho_wat : Water density, default is 1000 [kg/m³]
    rho_air : Air density, default is 1.2 [kg/m³]
    wat_visc : Water viscosity, default is 1e-3 [Pa s]


    """
    cfhalf = 2.25e-3
    if wind_speed < 0.1:
        cfhalf = 1.98e-3
    elif wind_speed < 3.06:
        cfhalf = 1.25e-3 * wind_speed**-0.2
    elif wind_speed < 22.3:
        cfhalf = (0.8+0.065*wind_speed)/1000

    friction_velocity = wind_speed * math.sqrt(rho_air/rho_wat) * math.sqrt(cfhalf)

    h = 0.01384 *(wind_speed*friction_velocity*rho_wat)/wat_visc

    b = 0

    if wind_speed < 5:
        b = 12.5 * schmdt_nmbr**0.667 + sct * math.log(schmdt_nmbr) / vnkarm - 5.3
    else:
        b = (0.55*math.sqrt(h)*(schmdt_nmbr**0.667 -0.2)- sct
             * math.log(schmdt_nmbr) / vnkarm + 11.2 * sct)

    height = (10 * friction_velocity * rho_wat) / wat_visc

    return 10*friction_velocity/(sct*math.log(height)/vnkarm+b+2.35)


def massic_flux_cohen(conc, part_coef, k,area):
    """
    Return the dissolution mass flux [kg/s]
    source : (Cohen et al., 1980)

    Parameters
    ----------
    conc : Component concentration in oil phase [kg/m³]
    part_coef : Partition coefficient ratio of component concentration
        (oil to water)[(g/cm³)/(g/cm³)]
    k : Water phase mass transfer coefficient (oil to water) [m/s]
    area : Area of the slick [m²]


    """
    return conc * part_coef*k*area



def massic_flux_shen(oil_sol, area, k):
    """
    Return the dissolution mass flux [kg/s]
    source : (Shen et al., 1993)

    Parameters
    ----------
    oil_sol : Oil solubility in water [kg/m³]
    area : Area of the slick [m²]
    k : Dissolution mass transfer coefficient [m/s]


    """
    return oil_sol * area * k



def molar_flux_HNS(k, area, sol, wat_conc,x = 1):
    """
    Return the dissolution molar flux [mol/s]
    source : (HNS-MS)

    Parameters
    ----------
    k : Mass transfert coefficient [m/s]
    area : Area of the slick [m²]
    sol : Pure component solubility [mol/m³]
    wat_conc : Water concentration of the component [mol/m³]
    x : molar fraction of the compounds, the default is 1


    """
    return k * area * (x * sol - wat_conc)
