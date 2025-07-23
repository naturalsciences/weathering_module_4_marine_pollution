# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 11:45:20 2021

@author: Ludovic Lepers
"""
import math



def interp(time, list_val, list_time):
    """
    Return the interpolation between the time and list_time for list_val.
    If outside of list_time, return the first or last element of list_val

    """
    if time <= list_time[0]:
        return list_val[0]
    elif time >= list_time[len(list_time)-1]:
        return list_val[len(list_val)-1]
    else:
        i = 0
        while list_time[i+1] < time:
            i += 1
    d1 = time - list_time[i]
    d2 = list_time[i+1]-time
    return (d2 * list_val[i]+d1 * list_val[i+1])/(d1+d2)



def stability(rho, mu, sat_cont, res_cont, asph_cont):
    """
    Return the stability
    source : (Fingas, 2015)

    Parameters
    ----------
    rho : Oil density [kg/m³]
    mu : Oil viscosity [Pa s]
    sat_cont : Saturate content [] (0.01-0.9)
    res_cont : Resin content [] (0.01-0.9)
    asph_cont : Asphaltene content [] (0.01-0.9)

    """
    den = math.exp(rho)
    visc = math.log(mu*1000)
    sst = 0
    if sat_cont < 0.45:
       sst = 45 -  sat_cont * 100

    else:
        sst = sat_cont * 100 - 45

    rst = 0
    if res_cont == 0:
        rst = 20
    elif res_cont < 10:
        rst = 10 - res_cont * 100
    else:
        rst = res_cont * 100 - 10


    ast = 0
    if asph_cont == 0:
        ast = 20
    elif asph_cont < 10:
        ast = 4 - asph_cont * 100
    else:
        ast = asph_cont - 4

    A_R = asph_cont / res_cont

    return (-15.3 + 1010*den - 3.66*visc + 0.174*rst - 0.579*ast + 34.4*A_R
            + 1.02*math.exp(den) - 7.91*math.exp(A_R) - 2740*math.log(den)
            + 12.2*math.log(visc) - 0.334*math.log(sst) - 3.17*math.log(rst)
            + 0.99*math.log(ast) - 2.29*math.log(A_R))





def emulsion_type(rho, visc, stability, asph, res):
    """
    Return the type of emulstion (entrained, unstable, stable, mesostable)

    Parameters
    ----------
    rho : Density [kg/m³]
    visc : Viscosity [Pa s]
    stability : Stability
    min_AR : Minimum between
    asph : Asphaltene content [] (0-1)
    res : Resine content [] (0-1)


    """
    min_AR = min(asph, res)
    if rho > 940 and visc > 6 and stability > -20 and stability < -3:
        return 'entrained'
    elif ((rho < 850 or rho > 1) and (visc < 0.1 or visc > 800)
          and (stability < -18 or stability > -4) and min_AR < 0.01):
        return 'unstable'
    elif stability > 4 and stability < 29:
        return 'stable'
    elif stability > -10 and stability < 5:
        return 'mesostable'
    else:
        raise Exception("Emulsion type not found")



def viscosity_from_emulsion_type(eml_type, time):
    """
    Return the viscosity coefficient from the emulsion stability []
    source : (Fingas, 2015)

    Parameters
    ----------
    eml_type : Emulsion type
    time : Time from the start of the emulsion [s]


    """
    time = time / (3600 * 24)
    list_visc = []
    if eml_type == 'entrained':
        list_visc = [1.9, 1.9, 2.1]
    elif eml_type == 'mesostable':
        list_visc = [7.2, 11, 32]
    elif eml_type == 'stable':
        list_visc = [405, 1054, 991]
    elif eml_type == 'unstable':
        list_visc = [0.99, 1, 1]

    return interp(time, list_visc, [1, 7, 365])



def water_cont_from_emulsion_type(eml_type, time):
    """
    Return the water cont from the emulsion stability []
    source : (Fingas, 2015)

    Parameters
    ----------
    eml_type : Emulsion type
    time : Time from the start of the emulsion [s]


    """
    time = time / (3600 * 24)
    list_visc = []
    if eml_type == 'entrained':
        list_visc = [0.445, 0.275, 0.06]
    elif eml_type == 'mesostable':
        list_visc = [0.643, 0.3, 0.06]
    elif eml_type == 'stable':
        list_visc = [0.81, 0.78, 0.7]
    elif eml_type == 'unstable':
        list_visc = [0.061, 0.06, 0.06]

    return interp(time, list_visc, [1, 7, 365])



def time_to_emulsion(eml_type, wave_height):
    """
    Return the time before the emulsion formation [s]
    source : (Fingas, 2015)

    Parameters
    ----------
    eml_type : Emulsion type
    wave_height : Wave height [m]


    """
    wave_height = wave_height *100
    param = []
    if eml_type == 'entrained':
        param = [27.1, 7.52]
    elif eml_type == 'mesostable':
        list_visc = [47, 49.1]
    elif eml_type == 'stable':
        list_visc = [30.8, 18.3]
    elif eml_type == 'unstable':
        raise Exception("Unstable stability do not lead to emulsion!")

    return 60*(param[0] + param[1]/wave_height**1.5)




def volume_OSERIT(surfc_amount, wave_height, k_em, dt, max_wat_cont = 0.8, c1 = 2000000):
    """
    Return the remaining surface volumic [m³/s] from OSERIT
    source : OSERIT manual

    Parameters
    ----------
    surfc_amount : Oil volume at the surface [m³]
    wave_height : Significant wave height [m]
    dt : timestep length [s]
    k_em : Oil ability to form emulsion (0-120)
        max_wat_cont : Maximum water content, the default is 0.8
    c1 : Constant, the default is 2000000



    """
    return surfc_amount*math.exp(-max_wat_cont/(1-max_wat_cont)*k_em/c1*wave_height*dt)






def wat_fract_MACKAY(wind_speed, fract_water, max_wat_cont = 0.8, c1 = 2e-6):
    """
    Return the emulsification fraction rate [1/s]
    source : (Fingas, 1995)

    Parameters
    ----------
    wind_speed : Wind speed [m/s]
    fract_water : Water content (0-1)
    max_wat_cont : Maximum water content, the default is 0.8
    c1 : COnstant, the default is 2e-5 or 2e-6 (contradictory sources)


    """
    return c1 * (wind_speed+1)**2*(1-fract_water/max_wat_cont)


def wat_fract_Eley(S_max, ks, water_dropplet_diam, time):
    """
    Return the water content of the oil []
    source: (Eley, 1988)

    Parameters
    ----------
    S_max : Maximum interfacial area [m²/cm³]
    ks : Interfacial parameter [1/s]
    water_dropplet_diam : Water dropplet diameter [m]
    time : Time from the start of the release [s]


    """
    S = S_max*(1-math.e**(ks*time))
    return S*water_dropplet_diam/(6+S*water_dropplet_diam)
