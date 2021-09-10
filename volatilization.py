# -*- coding: utf-8 -*-
"""
Created on Sep 07 2021

@author: Ludovic Lepers
"""
import math

def volatilization_coef(molar_weight, henry, temperature, R = 8.206e-5):
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


def henry(molar_weight, solubility, vapor_pressure):
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


def volatilization_rate(k, mass, dz, dt):
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
