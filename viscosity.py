# -*- coding: utf-8 -*-
"""
Created on Jan 26 09:42:32 2022

@author: Ludovic Lepers
"""

import math


def riazi_viscosity(A, B, I):
    """
    Return the viscosity [Pa s]
    source : (Riazi, 2001)

    params
    ------
    A : parameter [1/cp]
    B : parameter [1/cp]
    """

    return 1/(A + B /I) /1000

def riazi_A(Tbr, M, S, I20):
    """
    Return the parameter A [1/cp]
    source : (Riazi, 2001)

    params
    ------
    Tbr : Reduced boiling temperature []
    M : molar mass [kg/mol]
    S : specific gravity []
    I20 : refractive parameter at 20°C []
    """
    return (37.34745 - 0.20611 * M * 1000 + 141.1265 * S - 637.727 * I20
                - 6.75 * Tbr + 6.98 * Tbr * Tbr - 0.81 * Tbr * Tbr * Tbr)

def riazi_B(Tbr, M, S, I20):
    """
    Return the parameter B [1/cp]
    source : (Riazi, 2001)

    params
    ------
    Tbr : Reduced boiling temperature []
    M : molar mass [kg/mol]
    S : specific gravity []
    I20 : refractive parameter at 20°C []
    """
    return (-15.5437 + 0.046603 * M * 1000 - 42.8873 * S + 211.6542 * I20
                + 1.676 * Tbr - 1.8 * Tbr * Tbr + 0.212 * Tbr * Tbr * Tbr)

def riazi_reduced_Teb(Teb):
    """
    Return the reduced ebulition temperature []
    source : (Riazi, 2001)

    params
    ------
    Teb : Boiling temperature [K]
    """
    return (1.8 * Teb - 459.67) / 100

def riazi_I20(Teb, M, S):
    """
    Return the refraction parameter at 20°C [1/cp]
    source : (Riazi, 2001)

    params
    ------
    Teb : Boiling temperature [K]
    M : molar mass [kg/mol]
    S : specific gravity []
    """
    M = M * 1000 #to g/mol

    if M <= 300:
        return (2.3435e-2 * math.exp(7.029e-4 * Teb + 2.468 * S - 10.273e-4 * Teb * S)
                    * math.pow(Teb, 0.0572) * math.pow(S, -0.72))
    else:
        return (1.8429e-2 * math.exp(11.635e-4 * Teb + 5.144 * S - 5.92e-4 * Teb * S)
                    * math.pow(Teb,-0.4077) * math.pow(S, -3.333))


def riazi_I(T, I20):
    """
    Return the refraction parameter [] at the temperature T
    source : (Riazi, 2001)

    params
    ------
    T : temperature [K]
    I20 : refraction parameter [] at 20 °C
    """
    n20 = math.sqrt((2*I20+1)/(1-I20))
    n = n20 - 0.0004 *(T- 20 - 273.15)
    return (n * n - 1)/(n * n + 2)

def viscosity(T, Teb, M, rho):
    """
    Return the viscosity [Pa s]
    source : (Riazi, 2001)

    params
    ------
    T : Temperature [K]
    Teb : Boiling temperature [K]
    M : molar mass [kg/mol]
    rho : density [kg/m³]
    """
    S = rho / 1000
    I20 = riazi_I20(Teb, M, S)
    I = riazi_I(T, I20)
    Tbr = riazi_reduced_Teb(Teb)
    A = riazi_A(Tbr, M, S, I20)
    B = riazi_B(Tbr, M, S, I20)
    return riazi_viscosity(A, B, I)
