# -*- coding: utf-8 -*-
"""
Created on Sep 07 2021

@author: Ludovic Lepers
"""
import math


def henry(vapor_pressure, solubility, molar_mass, atm_pressure = 101325):
    """
    Return the henry constant [atm m³/mol]
    source : (Lyman et al., 1990)

    params
    ------
    vapor_pressure: vapor pressure [Pa]
    solubility: solubility [kg/m³]
    molar_mass: molar mass of the compounds [kg/mol]
    atm_pressure: atmospheric pressure, default if 101325 [Pa]
    """
    return (vapor_pressure / atm_pressure) / (solubility / molar_mass)


def nondimensional_henry(henry, temperature, R = 8.2e-5):
    """
    Return the nondimensional henry constant []
    source : (Lyman et al., 1990)

    params
    ------
    henry: henry constant [atm m³/mol]
    temperature: temperature [K]
    R: gaz constant 8.2e-5 [atm m³/mol]
    """
    return henry / (temperature * R)


def liquid_phase_exchange_coef(molar_mass, wind_speed):
    """
    Return the liquid phase exchange coefficent [cm/hr]
    source : (Lyman et al., 1990)

    params
    ------
    wind_speed: wind speed [m/s]
    molar_mass: molar mass [kg/mol]
    """
    molar_mass = molar_mass * 1000
    if(molar_mass < 65):
        return 20 * math.sqrt(44/molar_mass)
    elif wind_speed < 3:
        return 2.5
    elif wind_speed < 6:
        return 10
    elif wind_speed < 10:
        return 23
    else:
        return 50


def gas_phase_exchange_coef(molar_mass, wind_speed, current_speed):
    """
    Return the gas phase exchange coefficent [cm/hr]
    source : (Lyman et al., 1990)

    params
    ------
    wind_speed: wind speed [m/s]
    molar_mass: molar mass [kg/mol]
    current_speed: speed of the current [m/s]
    """
    molar_mass = molar_mass * 1000
    if(molar_mass < 65):
        return 3000 * math.sqrt(18/molar_mass)
    else:
        return 1137.5*(wind_speed+current_speed) * math.sqrt(18/molar_mass)


def mass_transfer_coefficient(nd_henry, gas_coef, liquid_coef):
    """
    Return the mass transfer coefficient [m/s]
    source : (Lyman et al., 1990)

    params
    ------
    nd_henry: wind speed []
    gas_coef: gas phase exchange coefficient [cm/hr]
    liquid_coef: liquid phase exchange coefficient [cm/hr]
    """
    a = 1 / (100 * 3600) #from cm/hr to m/s
    return a * (nd_henry * gas_coef * liquid_coef) / (nd_henry * gas_coef + liquid_coef)


def mass_flux_lyman(K, C, H, molar_mass, P = 0, atm_pressure = 101325):
    """
    Return the mass flux [kg/m² s]
    source : (Lyman et al., 1990)

    params
    ------
    K: mass transfer coefficent [m/s]
    C: concentration [kg/m³]
    H: henry constant [atm m³/mol]
    molar_mass: molar mass [kg/m³]
    P: vapor pressure (default is neglected at 0)[Pa]
    atm_pressure: atmospheric pressure, defaut is 101325 [Pa]
    """
    psh = P/(H * atm_pressure / molar_mass)
    return K * (C - psh)


def volatilization_coef_chemmap(molar_weight, henry, temperature, R = 8.206e-5):
    """
    Return volatilization coefficient [m/s]
    source : (“CHEMMAP technical User’s manual 6.10,” 2014)

    params
    ------
    molar_weight[kg/mol] (in the ref, the dimension is not defined)
    henry: henry constant
    temperature: temperature[K]
    R : gaz constant [atm m³ K / mole]
    """
    if henry < 3e-7:
        raise Exception("Volatilization should be neglected with this constant")
    h = henry /(R*temperature)

    mw = molar_weight / 1000
    a = 20 * math.sqrt(44/mw)
    b = 300 * math.sqrt(18/mw)
    k = (h*a*b)/(h*b+a) /(100*3600) #from cm/hr to m/s

    return k


def henry_chemmap(molar_weight, solubility, vapor_pressure):
    """
    Return the henry's law constant []. if < 3e-7, can be neglected
    source : (“CHEMMAP technical User’s manual 6.10,” 2014)

    params
    ------
    molar_weight[kg/mol] (in the ref, the dimension is not defined)
    solubility [kg/m³]
    vapor_pressure [Pa]
    """
    P = vapor_pressure/101325
    s = solubility * 1000
    mw = molar_weight / 1000

    return P / (s/mw)


def volatilization_rate_chemmap(k, mass, dz, dt):
    """
    Compute the transfer of the volatilization (kg/s)
    source : (“CHEMMAP technical User’s manual 6.10,” 2014)

    params
    ------
    k : volatilization coefficient [m/s]
    dz : vertical diffusivity [m²/s]
    dt : length of a timestep [s]
    mass: mass of the pollutant [kg]
    """

    return k * mass / math.sqrt(2 * dz * dt)
