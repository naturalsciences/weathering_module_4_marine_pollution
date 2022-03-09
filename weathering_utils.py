# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 13:41:10 2021

@author: Ludovic Lepers
"""
import evaporation as ev
import dissolution as di
import photooxidation as ph
import biodegradation as bi
import emulsification as em
import volatilization as vo
import oil_utils as otl
import matplotlib.pyplot as plt
import numpy as np
import copy
import math
<<<<<<< HEAD
import numpy as np

=======
import csv
>>>>>>> 47ab9a04cb38f60ea24a9c5e2eae49df8457476a

MAX_EVAPORATIVE_TEMP = 573.15


def to_half_life(days):
    """
    Return the constant [1/s] from the half life length [day]

    """
    s= 24 * 3600 * days
    return -math.log(1/2)/s

def create_csv(mix, matrix, name, path = ""):
    """
    Create a csv file with the results of the simulation

    params
    ------
    mix : mix represented by the matrix
    matrix : matrix ordered following (with the last row of comp being the sum):
        time [s]
        slick [m³]
        evaporated [m³]
        dissolved [m³]
        biodegraded [m³]
        photooxided [m³]
        emulsified [m³]
    name: name of the file (without .csv)
    path: path of the file if not in the same directory as the script
    """
    ind_tot = len(mix.list_component)
    with open(path+name+'.csv','w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writerow(['time [s]', 'remaining [m³]','evaporated [m³]','dissolved [m³]',
                            'biodegraded [m³]','photooxided [m³]','emulsified [m³]'])
        for i in range(len(matrix[:,0,0])):
            csv_writer.writerow([matrix[i,0,ind_tot],matrix[i,1,ind_tot],matrix[i,2,ind_tot],
                                matrix[i,3,ind_tot],matrix[i,4,ind_tot],matrix[i,5,ind_tot],
                                matrix[i,6,ind_tot]])



def plot_matrix_mix(mix, matrix, title_add = "", percent = False):
    """
    Plot the matrix

    Parameters
    ----------
    mix : mix represented by the matrix
    matrix : matrix ordered following (with the last row of comp being the sum):
        time [s]
        slick [m³]
        evaporated [m³]
        dissolved [m³]
        biodegraded [m³]
        photooxided [m³]
        emulsionned [m³]
    title_add: string which will be added to the title of the plot
    percent : if True, show the ordinate in percent, by default False

    """
    ind_tot = len(mix.list_component)
    sim_length = matrix[len(matrix)-1,0,0]
    sum = (matrix[0,1,ind_tot]+matrix[0,2,ind_tot]+matrix[0,3,ind_tot]
            +matrix[0,4,ind_tot]+matrix[0,5,ind_tot]+matrix[0,6,ind_tot])
    if percent:
        plt.plot(matrix[:,0,0]/3600,matrix[:,1,ind_tot]/sum*100)
        plt.plot(matrix[:,0,0]/3600,matrix[:,2,ind_tot]/sum*100)
        plt.plot(matrix[:,0,0]/3600,matrix[:,3,ind_tot]/sum*100)
        plt.plot(matrix[:,0,0]/3600,matrix[:,4,ind_tot]/sum*100)
        plt.plot(matrix[:,0,0]/3600,matrix[:,5,ind_tot]/sum*100)
        plt.plot(matrix[:,0,0]/3600,matrix[:,6,ind_tot]/sum*100)
        plt.ylabel('% (Volume)')
        plt.ylim([0, 100])
        #plt.ylim([0, matrix[0,1,ind_tot]]/sum*100)
    else :
        plt.plot(matrix[:,0,0]/3600,matrix[:,1,ind_tot])
        plt.plot(matrix[:,0,0]/3600,matrix[:,2,ind_tot])
        plt.plot(matrix[:,0,0]/3600,matrix[:,3,ind_tot])
        plt.plot(matrix[:,0,0]/3600,matrix[:,4,ind_tot])
        plt.plot(matrix[:,0,0]/3600,matrix[:,5,ind_tot])
        plt.plot(matrix[:,0,0]/3600,matrix[:,6,ind_tot])
        plt.ylabel('Volume [m³]')
        plt.ylim([0, matrix[0,1,ind_tot]])

    plt.xlim([0, sim_length/3600])
    plt.xlabel('Time[h]')
    plt.legend(['Remaining','Evaporated','Dissolved','Biodegradation','Photooxidation','Emulsion'])
    plt.title(title_add+mix.name)

    plt.show()



def compute_weathering(mix, temperature, wind_speed, sim_length, dt, water_volume,
                    slick_thickness, fix_area=None, apply_evaporation = 1,
                    apply_emulsion = 0, apply_volatilization = 0, apply_dissolution = 0,
                    wave_height = 0, vertical_diff = 1e-3, stability_class='C',
                    current_speed = 0, speed_for_dis=1e-3, start_dis = False):
    """
    Return a matrix with the amount for each timestep for each process.
    The mix is consider to be the remaining amount as slick. All the
    component with an ebulition point above or equals to 1000K will not
    evaporate.
    The matrix contains for each timestep and each pseudo component + the sum
    time [s] (same for all component)
    slick [m³]
    evaporated [m³]
    dissolved [m³]
    biodegraded [m³]
    photooxided [m³]
    emulsionned [m³]

    Parameters
    ----------
    mix : Mix to weather. Be carefull, it will be modified!
    temperature : Temperature [K]
    wind_speed : Wind speed [m/s]
    sim_length : Length of the simulation [s]
    dt : Length of an iteration [s]
    water_volume : Amount of water in which the mix is dissolved [m³]
    slick_thickness : Slick thickness [m], computed from area if it is not None
    fix_area : Vector with the size of the slick area [m²] at each timestep,
                if none, it will be computed from the slick_thickness and the volume.
                The default is None. Can be a scalar which is then transform into
                a vector
    stability_class : Pasquill stability index. The default is 'F'.
    apply_evaporation : 0 means no evaporation. If 1, use oiltrans for evaporation,
                    if 2 it will use ALOHA, if 3 it will use Fingas. The default is 1.
    apply_emulsion : 0 means no emulsion, 1 means emulsion as OSERIT, 2 as
                    Mackay. The default is 0.
    apply_volatilization : 0 means no volatilization, 1 means volatilization,
                    The default is 0.
    apply_dissolution : 0 means no dissolution, 1 means dissolution,
                    The default is 0.
    wave_height : Wave height [m]
    current_speed : Current speed [m/s]
    speed_for_dis: the speed of the molecule for the dissolution[m/s]
    start_dis: if True, will start dissolved and not in slick
    """
    #In case of an area which vary trought time
    if type(fix_area) != list and type(fix_area) != np.ndarray:
        if fix_area != None:
            area = fix_area
            fix_area = []
            for i in range(0, int(sim_length/dt)):
                fix_area.append(area)
    time_step_amount = int(sim_length / dt)
    matrix = np.zeros((int(time_step_amount),7, len(mix.list_component)+1))
    #allocation of the components in the matrix
    for i in range(0, len(mix.list_component)):
        if start_dis :
            matrix[0,3,i]=mix.get_comp(i).amount
            mix.get_comp(i).amount = 0
        else:
            matrix[0,1,i]=mix.get_comp(i).amount
    if start_dis :
        matrix[0,3,len(mix.list_component)]=sum(matrix[0,3])
    else:
        matrix[0,1,len(mix.list_component)]=sum(matrix[0,1])

    if apply_emulsion > 0:
        max_wat = mix.max_water
    for i in range(1,time_step_amount):
        matrix[i,0,:] = i * dt

        ev_fl = np.zeros(len(mix.list_component)+1)
        dis_fl = np.zeros(len(mix.list_component)+1)
        bio_fl = np.zeros(len(mix.list_component)+1)
        phot_fl = np.zeros(len(mix.list_component)+1)
        vol_fl = np.zeros(len(mix.list_component)+1)
        emul_fl = np.zeros(len(mix.list_component)+1)
        if fix_area is not None:
            area = fix_area[i]
            slick_thickness = mix.get_prop(temperature).amount / area
        else:
            area = mix.get_prop(temperature).amount / slick_thickness
        length = 2 * math.sqrt(area / math.pi)

        if mix.get_prop(temperature).amount> 0:

            schmdt_nmbr_air = ev.schmdt_nmbr_MW(mix.get_prop(temperature).molar_weight)#the one of chemmap can be used

            for j in range(0, len(mix.list_component)):
                comp = mix.get_comp(j)
                comp_area = comp.amount / matrix[i-1,1,len(mix.list_component)]* area

                if comp.amount > 0 : #for the volatilization
                    fract = mix.get_molar_fract(comp)
                    if comp.boiling_T < 1000 and ev.find_vapor_pressure(comp, temperature, MAX_EVAPORATIVE_TEMP) > 0:
                    #evaporation pseudo component
                        flux = 0
                        #OILTRANS
                        if apply_evaporation == 1 :
                            p_oil = ev.find_vapor_pressure(comp, temperature, MAX_EVAPORATIVE_TEMP)
                            molar_v = comp.molar_volume
                            #k from mw
                            k_evp = ev.mass_transfer_coefficient_OILTRANS(wind_speed, length,
                                                  schmdt_nmbr_air)
                            flux = ev.evap_volume_OILTRANS(k_evp, area, p_oil,
                                                          molar_v,temperature,
                                                          molar_fraction = fract)
                        #ALOHA
                        elif apply_evaporation == 2:
                            n = ev.pasquill_stability(stability_class)
                            wind_friction = ev.wind_friction_velocity(wind_speed, n)
                            Re = ev.roughness_reynolds(wind_friction)
                            kc = ev.graham_law(comp.molar_weight)
                            sc = ev.schmdt_nmbr_laminar(kc)

                            k = ev.mass_transfer_coefficient_ALOHA(length, Re, sc, n)


                            cs = ev.vapor_phase_sat_conc(comp.molar_weight,
                                                          ev.find_vapor_pressure(comp, temperature,
                                                          MAX_EVAPORATIVE_TEMP), temperature)

                            flux = (ev.evap_mass_flux_ALOHA(cs, wind_friction, k)
                                    * comp_area / comp.density)
                        flux = flux * dt
                        if comp.amount < flux:
                                flux = comp.amount
                        comp.amount -= flux
                        ev_fl[j] = flux

                    #dissolution
                        if apply_dissolution:
                            if (comp.solubility is not None
                                    and comp.molar_weight is not None):
                                Dc = di.diffusion_coefficient(comp.molar_volume)
                                schmdt_nmbr_water = ev.schmdt_nmbr(Dc)
                                #for surface

                                Sh = di.sherwood_slick(schmdt_nmbr_water, wind_speed, length)
                                k = di.mass_transfer_coefficient_HNS(Sh, Dc,length)

                                #for subsurface droplet
                                #d = 10**-5
                                #Sh = di.sherwood_droplet(schmdt_nmbr_water,0.001,d)
                                #k = di.mass_transfer_coefficient_HNS(Sh, Dc,d)

                                #k = di.mass_transfer_coefficient_HNS_Fern(wind_speed, schmdt_nmbr_water)
                                conc =  matrix[i-1,3,j]/(water_volume*comp.molar_volume)
                                flux = di.molar_flux_HNS(k, comp_area, comp.solubility / comp.molar_weight,
                                                           conc,fract)
                                flux = flux * comp.molar_volume *dt
                                if comp.amount < flux:
                                    flux = comp.amount

                                comp.amount -= flux
                                dis_fl[j] = flux

                    #biodegradation
                        if comp.h_l_biod is not None:
                            flux = bi.simple_half_live(matrix[i-1,3,j], comp.h_l_biod, dt)
                            if comp.amount < flux:
                                flux = comp.amount
                            comp.amount -= flux
                            bio_fl[j]=flux
                    #photooxidtion
                        if comp.h_l_phot is not None:
                            flux = ph.simple_half_live(comp.amount, comp.h_l_phot, dt)
                            if comp.amount < flux:
                                flux = comp.amount
                            comp.amount -= flux
                            phot_fl[j]=flux

            #evaporation fingas
            if apply_evaporation == 3:
                c1 = mix.fingas1
                c2 = mix.fingas2
                fract = ev.evap_fract_time_fingas(True, c1, matrix[i,0,0], temperature,c2) /100
                flux = fract * matrix[0,1,len(mix.list_component)] - matrix[i-1,2,len(mix.list_component)]
                if mix.get_prop(temperature).amount < flux:
                    flux = mix.get_prop(temperature).amount

                ev_fl[0:len(mix.list_component)] = -mix.add_amount(-flux)



            #emulsion: traited on all at once
            if apply_emulsion > 0 and mix.get_prop(temperature).amount > 0:
                #OSERIT
                if apply_emulsion == 1:
                    eml_rate = em.wat_volume_OSERIT(mix.get_prop(temperature).amount,wave_height, mix.K_em)*dt
                    flux = eml_rate/max_wat - eml_rate
                    if mix.get_prop(temperature).amount < flux:
                        flux = mix.get_prop(temperature).amount
                    emul_fl[0:len(mix.list_component)] = mix.add_amount(-flux)
                #MACKAY
                elif apply_emulsion == 2:
                    wat_amount = matrix[i-1,6,len(mix.list_component)] /(1-max_wat)*max_wat
                    wat_fract = wat_amount / (wat_amount+ mix.get_prop(temperature).amount +
                                              matrix[i-1,6,len(mix.list_component)])
                    eml_rate = em.wat_fract_MACKAY(wind_speed, wat_fract, max_wat)

                    eml_wat = eml_rate * dt
                    flux = (eml_wat/max_wat * (1-max_wat)
                            * ( matrix[i-1,6,len(mix.list_component)] +  mix.get_prop(temperature).amount))
                    if mix.get_prop(temperature).amount < flux:
                        flux = mix.get_prop(temperature).amount
                    emul_fl[0:len(mix.list_component)] = -mix.add_amount(-flux)

        #volatilization (after because comp.amount can be ==0)
        if apply_volatilization > 0:
            for j in range(0, len(mix.list_component)):
                comp = mix.get_comp(j)
                if matrix[i-1,3,j] > 0:
                    if (comp.get_partial_P(temperature) is not None
                            and comp.solubility is not None
                            and comp.molar_weight is not None):
                        flux = 0
                        if apply_volatilization == 1:
                            henry = vo.henry(comp.get_partial_P(temperature), comp.solubility,
                                                comp.molar_weight)
                            if henry > 3e-7:
                                nd_henry = vo.nondimensional_henry(henry, temperature)
                                kl = vo.liquid_phase_exchange_coef(comp.molar_weight, wind_speed)
                                kg = vo.gas_phase_exchange_coef(comp.molar_weight, wind_speed, current_speed)
                                K = vo.mass_transfer_coefficient(nd_henry, kg, kl)
                                conc =  matrix[i-1,3,j] * comp.density/(water_volume)
                                flux = vo.mass_flux_lyman(K, conc, henry, comp.molar_weight)
                                flux = flux * dt * area / comp.get_density(temperature)

                        elif apply_volatilization == 2:
                                henry = vo.henry_chemmap(comp.molar_weight, comp.solubility,
                                                        comp.get_partial_P(temperature))

                                if henry >= 3e-7:
                                    k = vo.volatilization_coef_chemmap(comp.molar_weight,
                                                                henry, temperature)
                                    rate = vo.volatilization_rate_chemmap(k, matrix[i-1,3,j], 1e-3, dt)
                                    flux = rate*dt / comp.get_density(temperature)

                        if flux > matrix[i-1,3,j]:
                            flux = matrix[i-1,3,j]
                        vol_fl[j] = flux

        matrix[i,1,0:len(mix.list_component)] = mix.get_array_amount()
        matrix[i,2] = ev_fl + matrix[i-1,2] + vol_fl
        matrix[i,3] = dis_fl  + matrix[i-1,3] - vol_fl
        matrix[i,4] = bio_fl  + matrix[i-1,4]
        matrix[i,5] = phot_fl  + matrix[i-1,5]
        matrix[i,6] = emul_fl  + matrix[i-1,6]
        #sums

        for j in range(1,7):
            matrix[i,j,len(mix.list_component)] = sum(matrix[i,j,0:len(mix.list_component)])
    return matrix



def return_oils(amount, temperature):
    """
    Returns a list of some oils

    Parameters
    ----------
    amount : Amount of each mix [m³]
    temperature : Temperature [K]

    """
    oils_pseudo_comp = []

    brent_blend = otl.mix('BRENT BLEND')
    brent_blend_cut_T = [40,80,100,120,150,160,180,200,250,300,350,400,500,600,700]
    brent_blend_fract = [3.0,4.0,5.0,19.0,22.0,25.0,29.0,32.0,42.0,52.0,62.0,70.0,85.0,95.0,99.0]
    brent_blend.generate_component_cut(brent_blend_cut_T, brent_blend_fract, amount)
    brent_blend.density = 835
    brent_blend.K_em = 20
    brent_blend.visco = 4.5
    brent_blend.add_Fingas(3.39, 0.048)
    oils_pseudo_comp.append(brent_blend)


    arabian = otl.mix('ARABIAN')
    arabian_cut_T = [69,100,128,151,172,198,213,238,261,282,304,329,351,373,398]
    arabian_fract = [5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0,75.0]
    arabian.density = 863
    arabian.add_Fingas(1.89)
    arabian.K_em = 100
    arabian.visco = 40
    arabian.generate_component_cut(arabian_cut_T, arabian_fract, amount)
    oils_pseudo_comp.append(arabian)



    alba = otl.mix('ALBA')
    alba_cut_T = [20,149,204,260,277,290,303,311,322,333,343,385,454,485,534]
    alba_fract = [1.0,4.0,9.0,11.0,14.0,17.0,19.0,21.0,23.0,26.0,28.0,36.0,50.0,55.0,67.0]
    alba.density = 934
    alba.add_Fingas(0.450)
    alba.K_em = 20
    alba.visco = 463.5
    alba.generate_component_cut(alba_cut_T, alba_fract, amount)
    oils_pseudo_comp.append(alba)


    boscan = otl.mix('BOSCAN')
    boscan_cut_T = [91,230,287,350,369,384]
    boscan_fract = [1.0,5.0,10.0,20.0,30.0,40.0]
    boscan.density = 999
    boscan.add_Fingas(-0.15, 0.013)
    boscan.K_em = 1e-4
    boscan.visco = 3
    boscan.generate_component_cut(boscan_cut_T, boscan_fract, amount)
    oils_pseudo_comp.append(boscan)

    for oil in oils_pseudo_comp:
        add_oil_properties(oil, temperature)

    return oils_pseudo_comp


def return_chemicals(amount):
    """
    Return a list of some chemicals


    """
    chemicals_list = []

    butyl_acetate = otl.mix('Butyl acetate')
    ba = otl.component('Bulk',amount)#HNS MS 20°C
    ba.density = 881
    ba.molar_weight = 0.1162
    ba.boiling_T = 126+273.15
    ba.partial_P = 1990
    ba.ref_T_Clau = 25+273
    ba.vap_enthalpie = 377582
    ba.solubility = 5670e-3
    ba.h_l_biod = to_half_life(10)
    butyl_acetate.add_component(ba)

    chemicals_list.append(butyl_acetate)



    butyl_acrylate = otl.mix('Butyl acrylate')
    bay = otl.component('Bulk',amount)#HNS MS 20°C
    bay.density = 900
    bay.molar_weight = 0.12817
    bay.boiling_T = 148.8+273.15
    bay.partial_P = 727
    bay.ref_T_Clau = 25+273
    bay.vap_enthalpie = 458981 #by computing from the partial pressure at 20 and 25
    bay.solubility = 1440e-3
    bay.h_l_biod = to_half_life(10)
    bay.h_l_phot = to_half_life(0.0138)
    butyl_acrylate.add_component(bay)

    chemicals_list.append(butyl_acrylate)


    #2ethylhexyl acrylate
    ethylhexyl_acrylate = otl.mix('2-ethylhexyl acrylate')
    eay = otl.component('Bulk',amount)#HNS MS 20°C
    eay.density = 890
    eay.molar_weight = 0.18428
    eay.boiling_T = 216+273.15
    eay.partial_P = 24
    eay.ref_T_Clau = 25+273
    eay.vap_enthalpie = 397003
    eay.solubility = 21e-3
    eay.h_l_biod = to_half_life(10)
    ethylhexyl_acrylate.add_component(eay)

    chemicals_list.append(ethylhexyl_acrylate)


    heptane = otl.mix('Heptane')
    hep = otl.component('Bulk',amount)#HNS MS 20°C
    hep.density = 680
    hep.molar_weight = 0.1002
    hep.boiling_T = 98+273.15
    hep.partial_P = 6133
    hep.ref_T_Clau = 25+273
    hep.vap_enthalpie = 364960
    hep.solubility = 0.48e-3
    hep.h_l_biod = to_half_life(10)
    heptane.add_component(hep)

    chemicals_list.append(heptane)


    toluene = otl.mix('Toluene')
    tol = otl.component('Bulk',amount)#HNS MS 20°C
    tol.density = 868.3
    tol.molar_weight = 0.09215
    tol.boiling_T = 110.58+273.15
    tol.partial_P = 3800
    tol.ref_T_Clau = 25+273
    tol.vap_enthalpie = 412480
    tol.solubility = 110e-3
    tol.h_l_biod = to_half_life(30)

    tol.mu_max = 1.5e-4
    tol.ks = 1.96e-3
    tol.Y_oil = 1.22


    toluene.add_component(tol)

    chemicals_list.append(toluene)


    methanol = otl.mix('Methanol')
    meth = otl.component('Bulk',amount)#HNS MS 20°C
    meth.density = 791.4
    meth.molar_weight = 0.032042
    meth.boiling_T = 64.6+273.15
    meth.partial_P = 12265
    meth.ref_T_Clau = 20+273
    meth.vap_enthalpie = 1168154
    meth.solubility = 10000 #at 0% of salt
    methanol.add_component(meth)

    chemicals_list.append(methanol)



    pentane = otl.mix('Pentane')
    pen = otl.component('Bulk',amount)#HNS MS 20°C
    pen.density = 626.2
    pen.molar_weight = 0.072149
    pen.boiling_T = 36.06+273.15
    pen.partial_P = 57328
    pen.ref_T_Clau = 20+273
    pen.vap_enthalpie = 366190
    pen.solubility = 27e-3
    pentane.add_component(pen)

    chemicals_list.append(pentane)



    acetone = otl.mix('Acetone')
    ace = otl.component('Bulk',amount)#HNS MS 20°C
    ace.density = 790
    ace.molar_weight = 0.05808
    ace.boiling_T = 56.2+273.15
    ace.partial_P = 30930
    ace.ref_T_Clau = 25+273
    ace.vap_enthalpie = 533574
    ace.solubility = 790 #at 0% of salt
    tol.h_l_biod = to_half_life(10)
    tol.h_l_phot = to_half_life(1.67)
    acetone.add_component(ace)

    chemicals_list.append(acetone)



    xylene = otl.mix('Xylene')
    xyl = otl.component('Bulk',amount)#HNS MS 20°C
    xyl.density = 870
    xyl.molar_weight = 0.10616
    xyl.boiling_T = 140.2+273.15
    xyl.partial_P = 1070
    xyl.ref_T_Clau = 25+273
    xyl.vap_enthalpie = 401733
    xyl.solubility = 100e-3
    xylene.add_component(xyl)

    chemicals_list.append(xylene)

    return chemicals_list
