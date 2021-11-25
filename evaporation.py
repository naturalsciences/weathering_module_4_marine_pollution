# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 08:49:11 2021

@author: Ludovic Lepers
"""

import numpy as np
import math


def density_to_API(rho15):
    """
    Return the API from the density

    Parameters
    ----------
    rho15 : Density at 15°C [kg/m³]

    """
    d = rho15 / 999.06      #cf https://www.engineeringtoolbox.com/water-density-specific-weight-d_595.html?vA=15&units=C#
    return (141.5/d) - 131.5



def API_to_density(API):
    """
    Return the density from the API

    Parameters
    ----------
    API

    """
    if API < 0:
        API = 0
    d = 141.5/(API+131.5)
    return d * 999.06  #cf https://www.engineeringtoolbox.com/water-density-specific-weight-d_595.html?vA=15&units=C#


def grad_dist_T_rho(rho15):
    """
    Return the gradient distililation temperature [K] from the density
    source : (Mishra and Kumar, 2015)

    Parameters
    ----------
    rho15 : Density at 15°C [kg/m³]

    """
    return 985.62 - 13.597 * density_to_API(rho15)



def boiling_T_rho(rho15):
    """
    Return the boiling temperature [K] from the density
    source : (Mishra and Kumar, 2015)

    Parameters
    ----------
    rho15 : Density at 15°C [kg/m³]

    """
    return 532.98 - 3.125 * density_to_API(rho15)


def rho_from_boiling_T(T):
    """
    Return the density [Kg/m³] from the boiling temperature [K]
    source : (Mishra and Kumar, 2015)

    Parameters
    ----------
    T : Boiling temperature [K]

    """
    return API_to_density((532.98-T) / 3.125)



def graham_law(MW1, MW2 = 0.018, diff2 = 2.39*10**-5):
    """
    Return the molecular diffusivity with the graham laws in air

    Parameters
    ----------
    MW1 : Molecular weight of the chemical from interest [kg/mol]
    MW2 : Molecular weight of the reference, here water and 18 [kg/mol]
    diff2 : Molecular diffusivity of the reference, here water at 2.39*10**-5 [m²/s]

    """
    if MW1 is None:
        raise Exception("The component have no molecular weight define")
    return diff2 * math.sqrt(MW2/MW1)



def schmdt_nmbr(diffusivity, kinematic_viscosity = 10**-6):
    """
    Return the Schmidt number []

    Parameters
    ----------
    kinematic_viscosity : kinematic viscosity of water (defautl 10^-6)[m²/s]
    diffusivity : mass diffusivity of the chemical in water [m²/s]


    """
    return kinematic_viscosity / diffusivity



def schmdt_nmbr_MW(MW):
    """
    Return the Schmidt number [] from the mole weighted average of the oil,
    usefull in case of multi-component
    source : (Berry et al.2012)

    Parameters
    ----------
    MW : Mole weighted average of the oil [kg/mol]


    """

    return 1.3676 * (0.018/(MW*1000))**(-1/2)



def schmdt_nmbr_MW_chemmap(MW):
    """
    Return the Schmidt number [] in air from the molar weight of the compounds,
    usefull in case of multi-component (from chemap)
    source : (“CHEMMAP technical User’s manual 6.10,” 2014)

    Parameters
    ----------
    MW : Mole weighted average of the oil [kg/mol]


    """
    Dref = 0.0556 #cumene
    MWref = 120
    if MW < 0.1: #for pentane
        Dref = 0.071
        MWref = 72.15


    return 0.15/(Dref * (MWref/(MW*1000))**-(1/2))


def schmdt_nmbr_laminar(kc, kinematic_visc = 1.516e-5):
    """
    Return the laminar Schmidt number [], it uses air coefficient as default
    source : ALOHA

    Parameters
    ----------
    kc : molecular diffusivity [m²/s]
    kinematic_visc : kinematic viscosity of air, default 1.516e-5 [m²/s] (20°C)

    """
    return kinematic_visc/kc


def mass_transfer_coefficient_george(D, u, L, MmS, T, p_oil,
                                     molar_volume, R = 8.314, vis = 1.48*10**-5,
                                     Mminf = 0.02896, atm_p = 101325 ):
    """
    Return the mass transfert coefficient alpha
    source : (Kotzakoulakis and George, 2018)

    Parameters
    ----------
    D : Diffusion coefficient in air [m²/s]
    u : Air velocity [m/s]
    L : Length of the spill in the wind direction [m]
    MmS : Molar mass of air and vapor right above the spill [kg/mol]
    T : Temperature [K]
    p_oil : Vapor pressure of the oil [Pa]
    molar_volume : Molar volume of the liquid at the surface [m³/mol]
    R : Perfect gas constant, the default is 8.314 [J/mol K]
    vis : Kinematic viscosity of air : 1.48*10**-5 [m²/s]
    Mminf : Molar mass of air away from the spill, default is 0.02896 [kg/mo]
    atm_P : Atmospheric pressure [Pa] the default is 101 325 Pa


    """
    log_p = p_oil / math.ln(atm_p/(atm_p-p_oil))
    a = MmS / Mminf
    b = -0.662 * D**(2/3) * vis**(-1/6) * math.sqrt(u * L) / L * math.ln(a)
    return b / (a-1) * atm_p / log_p * p_oil * molar_volume / (R * T)



def mass_transfer_coefficient_mackay_matsugu(wind_speed, diameter, schmdt_nmbr,
                                             n = 0.25):
    """
    Return the mass transfer coefficient [m/s]
    source : (Mackay and Matsugu, 1973)


    Parameters
    ----------
    wind_speed : Wind speed 10 meters above the surface [m/s]
    diameter : Pool diameter [m]
    schmdt_nmbr : Schmidt number []
    n : Wind profile (0.25->1), the default is 0.25.

    """

    C = 0.0292 * schmdt_nmbr ** -0.47
    q = (2-n)/(2+n)
    r = -n/(2+n)

    return C * ((wind_speed * 3600) ** q) * (diameter ** r) / 3600


def mass_transfer_coefficient_mishra_kumar(wind_speed):
    """
    Return the mass transfert coefficient [m/s] only from wind speed
    source:(Berry et al., 2012)

    Parameters
    ----------
    wind_speed : Wind speed 10 meters above the surface [m/s]

    """
    return 0.0025 * (wind_speed**0.78)


def mass_transfer_coefficient_OILTRANS(wind_speed, length, schmdt_nmbr):
    """
    Return the mass transfert coefficien [m/s], formula derived
    from (Mackay and Matsugu, 1973)
    source : (Berry et al., 2012)

    Parameters
    ----------
    wind_speed : Wind speed 10 meters above the surface [m/s]
    length : Downwind length of the oil slick axis [m]
    schmdt_nmbr : Schmidt number []


    """
    return (0.0048 * (wind_speed **(7/9)) * (length ** (-1/9))
            * schmdt_nmbr ** (-2/3))



def mass_transfer_coefficient_ALOHA(Dp, Re, sc, n, k = 0.4, y = 0.577, SCt = 0.85,
                                    z = 0.0002):
    """
    Return the mass transfert coefficient from ALOHA []
    source : ALOHA

    Parameters
    ----------
    Dp : Diameter along the wind axis [m]
    Re : Roughness Reynolds number []
    sc : Laminar schmidt number in air []
    n :
    k : Von Karmon constant, equals to 0.4
    y : Euler's constant, equals to 0.577
    SCt : Turbulent schmidt number, equals to 0.85
    z : Surface roughness lenght [m]. The default for open sea is 0.0002 m (wikipedia)


    """

    X = (n*k**2*Dp)/(SCt*z*math.e**(1/n))

    smooth = (3.85*sc**(1/3)-1.3)**2 + SCt/k*math.log(0.13*sc)
    rough = 7.3*Re**(1/4)*math.sqrt(sc)-5*SCt

    fsc = 0

    if Re < 0.13:
        fsc = smooth
    elif Re > 2:
        fsc = rough
    else:
        fsc = smooth * (1-Re)+rough*Re

    lmbd = 1/n +1 + 2*math.log(1+n)-2*y+k/SCt*(1+n)*fsc


    pi = math.pi
    p1 = k/SCt*(1+n)
    lne = math.log(X*math.e**lmbd)
    p2 = 1/2-1/pi*math.atan(lne/pi)
    p3 = (1-y)/(lne**2+pi**2)
    p4 = ((1+(1-y)**2+1/6*pi**2)*lne)/(lne**2+pi**2)**2
    return p1*(p2+p3+p4)


def roughness_reynolds(wind_friction, kinematic_visc = 1.516e-5,z = 0.0002):
    """
    Return the roughness Reynolds number []
    source : ALOHA

    Parameters
    ----------
    wind_friction : Wind friction velocity [m/s]
    kinematic_visc :kinematic viscosity of air, default 1.516e-5 [m²/s] (20°C)
    z : Surface roughness lenght [m]. The default for open sea is 0.0002 m (wikipedia)

    """
    return (wind_friction*z)/kinematic_visc


def henry(p_oil, mol_vol, T, R = 8.314):
    """
    Return Henry constant []
    source : (Stiver and Mackay, 1984)

    Parameters
    ----------
    p_oil : Vapor pressure of the hydrocarbon [Pa]
    mol_vol : Liquid molar volume [m³/mol]
    T : Temperature [K]
    R : Perfect gas constant, the default is 8314 [J/mol K]


    """
    return (p_oil * mol_vol)/(T * R)


def henry_T(ebulition_T, T):
    """
    Return Henry constant [] from ebulition temperature
    source : (Mishra and Kumar, 2015)

    Parameters
    ----------
    ebulition_T : Ebulition temperature [K] (overpredicted above 280°C)
    T : Temperature [K]


    """
    if ebulition_T > (280+273.15):
        print("Ebulition temperature for computing henry above 280 °C !")

    return math.exp(6.3-10.3*(ebulition_T/T))


def theta_wtt(k, area, initial_volume):
    """
    Return the values of theta without the time [1/s]
    source: (Stiver and Mackay, 1984)

    Parameters
    ----------
    k : Mass transsfert coefficent [m/s]
    area : Area [m²]
    initial_volume : Initial spill volume [m³]

    """
    return (k * area) / initial_volume



def C_fingas(heavy, precent_dist_180):
    """
    Returns the constant C used in fingas model
    source : (Fingas, 2015)

    Parameters
    ----------
    heavy : True if the fuel need to follow a ln, else it will be a sqrt.
    precent_dist_180 : Percentage distilled by weight at 180 °C []


    """
    if heavy:
        return 0.165 * precent_dist_180
    else:
        return 0.0254 * precent_dist_180



def C_fingas_corr(C, thickness):
    """
    Returns the constant C used in fingas model corrected with the
    slick thickness, used if the slick is less than 1.5 mm
    source : (Fingas, 2015)

    Parameters
    ----------
    C : Fingas constant
    thickness : Slick thickness  [m]


    """
    if thickness > 1.5 / 1000:
        print("Slick thickness is more than 1.5 mm, correction not needed")

    return C + 1 - 0.78 * math.sqrt(thickness/1000)



def molar_volume_eb_T(eb_T):
    """
    Return the molar volume [m³/mol] from the ebulition temperature
    source :(Berry et al., 2012)

    Parameters
    ----------
    eb_T : Boiling point cut [K]


    """
    return 7 * 10**-5 - (2.102 * 10**-7 * eb_T) + 10**-9 * eb_T**2



def vapor_pressure_eb_T(eb_T, T, atm_P = 101325, R = 8.134):
    """
    Return the vapor pressure [Pa] from the boiling point cut
    source :(Berry et al., 2012)

    Parameters
    ----------
    eb_T : Boiling point cut [K]
    T : Temperature [K]
    atm_P : Atmospheric pressure [Pa] the default is 101 325 Pa
    R : Perfect gas constant, the default is 8.314 [J/mol K]

    """
    dsi = 8.75 + 1.987*math.log10(eb_T)
    c2 = 0.19 * (eb_T-18)

    a = (dsi * (eb_T-c2)**2 )/(R*eb_T)
    b = (1/(eb_T-c2)-1/(T-c2))
    print(eb_T,dsi, c2, a, b)
    return atm_P * math.exp(a*b)


def find_vapor_pressure(comp, temperature, max_ev_temp):
    """
    Compute the vapor pressure from the temperature, if not already computed
    Everything with a ebulition temperature above max_ev_temp does not evaporate
    if his vapor pressure is not already defined
    Parameters
    ----------
    comp : component
    temperature : Temperature [K]
    max_ev_temp : temperature [K] above which the partial pressure is put at 0

    """
    if comp.get_partial_P(temperature) is None:
        boiling_T = comp.boiling_T
        if boiling_T is None:
            boiling_T = boiling_T_rho(comp.density)
        p_oil = vapor_pressure_eb_T(boiling_T,temperature)
        if(boiling_T > max_ev_temp):
            p_oil = 0
    else:
        p_oil = comp.get_partial_P(temperature)

    return p_oil


def vapor_phase_sat_conc(MW, vapour_pressure, temperature, R=8.314):
    """
    Return the vapor phase saturation concentration [kg/m³]

    Parameters
    ----------
    MW : Molecular weight [kg/mol]
    vapour_pressure : Vapour pressure [Pa]
    T : Temperature [K]
    R : Perfect gas constant, the default is 8.314 [J/mol K]
    source:ALOHA

    """
    return MW * vapour_pressure /(temperature * R)


def wind_friction_velocity(wind_speed, n, height = 10):
    """
    Return the friction velocity [m/s]
    source : ALOHA

    Parameters
    ----------
    wind_speed : Speed of the wind [m/s] at the height
    n : Valaue of the exponent for Pasquill stability classes
    height : Height of wind measure, the default is 10 [m]


    """
    return wind_speed * 0.03 * (10/height)**n



def pasquill_stability(stability_class):
    """
    return the pasquill stability index
    source : ALOHA

    Parameters
    ----------
    stability_class : Stability class (A->F)

    """
    switcher = {
        'A': 0.108,
        'B': 0.112,
        'C': 0.120,
        'D': 0.142,
        'E': 0.203,
        'F': 0.253
    }
    return switcher.get(stability_class, None)



def massic_rate_sutton(k, conc_evaporating_fluid,
                       wind_speed, pool_area, schmdt_nmbr, r):
    """
    Return the evaporation rate [kg/s m²] of sutton's model
    source : (Fingas, 2015)

    Parameters
    ----------
    k : mass transsfert coefficent [1/m²]
    conc_evaporating_fluid : Concentration of the evaporating fluid [kg/m³]
    wind_speed : Wind speed 10 meters above the surface [m/s]
    pool_area : Area of the pool []
    schmdt_nmbr :  Schmidt number []
    r : empirical (0->2/3) []

    """
    return (k * conc_evaporating_fluid
         * (wind_speed**(7/9)) * (pool_area**(1/9)) * (schmdt_nmbr**r))


def molar_rate_mackay_matsugu(k, p_oil, T, R = 8.314, p_bulk = 0):
    """
    Return the molar evaporation rate [mol/m² s]
    source : (Mackay and Matsugu, 1973)

    Parameters
    ----------
    k : mass transsfert coefficent [m/s]
    p_oil : Vapor pressure of the hydrocarbon [Pa]
    T : Temperature [K]
    R : Perfect gas constant, the default is 8.314 [J/mol K]
    p_bulk : Vapor pressure of the hydrocarbon far from surface,default is 0 [Pa]


    """
    return k *((p_oil-p_bulk)/(R*T))



def evap_fract_stiver_mackay_integr(henry, theta_wtt):
    """
    Return the evaporated fraction by unit of time  [1/s],
    consider the area constant
    source:(Stiver and Mackay, 1984)

    Parameters
    ----------
    henry : Henry constant []
    theta_wtt : Theta without the time [1/s]

    """
    return henry * theta_wtt



def evap_fract_time_stiver_mackay(T, theta_wtt, time, k1, k2, k3):
    """
    Return the evaporated fraction []
    source:(Stiver and Mackay, 1984)

    Parameters
    ----------
    T : Temperature [K]
    theta_wtt : Theta without the time [1/s]
    time : Time from the spill [s]
    k1 : empirical constant
    k2 : empirical constant
    k3 : empirical constant

    """
    return T/k1 * math.log(1+(k1*theta_wtt*time)/T) * math.exp(k2 - k3/T)



def evap_fract_mishra_kumar(theta_wtt, evaporated_fract, T, eb_T, dist_T):
    """
    Return the evaporated fraction by unit of time [/s]
    source:(Stiver and Mackay, 1984)

    Parameters
    ----------
    theta_wtt : Theta without the time [1/s]
    evaporated_fract : Fraction of oil evaporated []
    T : Temperature [K]
    eb_T : Initial boiling temperature [K]
    dist_T : Gradient of oil distribution curve [K]


    """
    if dist_T < 0:
        print("Distilation temperature negative!")
    return (theta_wtt *
            math.exp(6.3-10.3*((eb_T + dist_T * evaporated_fract) /T)))



def evap_fract_time_fingas(heavy, c1, time, T, c2 = None):
    """
    Return the evaporated fraction []
    source : (Fingas, 2015)

    Parameters
    ----------
    heavy : True if the fuel need to follow a ln, else it will be a sqrt
    C : fingas constant
    time : Time from the spill [s]
    T : Temperature [K]
    c2 : second fingas constant, can be compute by defaut []

    """
    if c1 is None:
        raise Exception("The mix does not have Fingas constant c1 defined")

    if(time < 60): #else it will cause trouble with the log
        return 0
    b = 0
    T = T - 273.15
    if heavy:
        if c2 is None:
            c2 = 0.045
        #b = 15+273.15
        return (c1 + c2 * (T - b)) * math.log(time/60)  #temperature in [K], time in sec
    else:
        if c2 is None:
            c2 = 0.01
        #b = 15+273.15
        return (c1 + c2 * (T - b)) * math.sqrt(time/60)  #temperature in [K], time in sec



def evap_volume_OILTRANS(k, area, p_oil, molar_v, T, molar_fraction = 1,
                         R = 8.314):
    """
    Return the evaporated fraction by unit of time [/s]
    source : (Berry et al., 2012)

    Parameters
    ----------
    k : mass transsfert coefficent [m/s]
    area : Area [m²]
    p_oil : Vapor pressure of the hydrocarbon [Pa]
    molar_v : Molar volume [m³/mol]
    T : Temperature [K]
    molar_fraction : Molar fraction, can be other than 1 if not pure
    R : Perfect gas constant, the default is 8.314 [J/mol K]


    """
    return (k * area * p_oil * molar_v * molar_fraction) / (R * T)


def evap_mass_flux_ALOHA(cs, wind_friction, coef):
    """
    Return the evaporative mass flux [kg/m²s]
    source : ALOHA

    Parameters
    ----------
    cs : Chemical vapor phase saturation concentration [kg/m³]
    frict_vel : Friction velocity of the air [m/s]
    coef : Mass transfert coefficient []

    """
    return cs * wind_friction * coef
