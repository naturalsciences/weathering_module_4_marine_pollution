# -*- coding: utf-8 -*-

import evaporation as ev
import oil_utils as ou
import weathering_utils as wu  #contains the usefull fonctions
import iojson
import csv
import copy
import numpy as np
import random
from matplotlib import pyplot as plt


def createMixImaros_random(name, volume_init, list_proportion, density,t_density, viscosity, t_viscosity):
    """
    sum of list_proportion (it is the massic proportion)
    !!!!viscosity is in mPas
    """
    mix = ou.mix(name)
    mix.density = density
    mix.density_T = t_density
    if len(viscosity) > 0:
        for visc in viscosity:
             visc = visc/1000
        mix.viscosity = viscosity
        mix.viscosity_T = t_viscosity

    c1c4 = ou.component('c1-c4',list_proportion[0])
    #c1c4.density = 1.6433#gas
    c1c4.density = rand(504,573)#liquid
    c1c4.molar_weight = rand(0.016043,0.058)
    c1c4.boiling_T = rand(112.15,273.15)
    c1c4.partial_P = rand(170e3,3.84e6)
    c1c4.solubility = rand(22e-3,61e-3)
    c1c4.freezing_T = rand(85.5,136)
    c1c4.compute_molar_volume()
    mix.add_component(c1c4)

    c5 = ou.component('c5 sat',list_proportion[1])
    c5.density = rand(616,751)
    c5.molar_weight = rand(0.07,0.072)
    c5.boiling_T = rand(301,322.3)
    c5.partial_P = rand(45000,76992)
    c5.ref_T_Clau = rand(293.15,293.15)
    c5.vap_enthalpie = rand(366190,366190)
    c5.solubility = rand(27e-3,156e-3)
    c5.freezing_T = rand(113,179.2)
    c5.compute_molar_volume()
    mix.add_component(c5)

    c6 = ou.component('c6 sat',list_proportion[2])
    c6.density = rand(653,778.1)
    c6.molar_weight = rand(0.084159,0.08618)
    c6.boiling_T = rand(334,353.88)
    c6.partial_P = rand(15999,19998)
    c6.ref_T_Clau = rand(298.15,303.15)
    c6.vap_enthalpie = rand(366295,390810)
    c6.solubility = rand(6.6e-3,45e-3)
    c6.freezing_T = rand(120, 279.74)
    c6.compute_molar_volume()
    mix.add_component(c6)

    c7 = ou.component('c7 sat',list_proportion[3])
    c7.density = rand(679,693.8)
    c7.molar_weight = rand(0.098,0.1)
    c7.boiling_T = rand(362.85,374)
    c7.partial_P = rand(6133,6133)
    c7.ref_T_Clau = rand(298.15,298.15)
    c7.vap_enthalpie = rand(364960,364960)
    c7.solubility = rand(0.48e-3,0.014)
    c7.freezing_T = rand(146.8, 182.15)
    c7.compute_molar_volume()
    mix.add_component(c7)

    c8 = ou.component('c8 sat',list_proportion[4])
    c8.density = rand(696,698)
    c8.molar_weight = rand(0.114,0.114)
    c8.boiling_T = rand(391,399.15)
    c8.partial_P = 5300
    c8.ref_T_Clau = 310.15
    c8.solubility = 0.66e-3
    c8.freezing_T = 152
    c8.compute_enthalpy()
    c8.compute_molar_volume()
    mix.add_component(c8)

    c9 = ou.component('c9 sat',list_proportion[5])
    c9.density = 718
    c9.molar_weight = 0.128
    c9.boiling_T = 424
    c9.partial_P = 590
    c9.ref_T_Clau = 298.15
    c9.solubility = 0 #put to 0
    c9.freezing_T = 219.5
    c9.compute_enthalpy()
    c9.compute_molar_volume()
    mix.add_component(c9)

    benzene = ou.component('benzene',list_proportion[6])
    benzene.density = 880
    benzene.molar_weight = 0.07811
    benzene.boiling_T = 353.25
    benzene.partial_P = 12700
    benzene.ref_T_Clau = 298.15
    benzene.solubility = 1.360
    benzene.vap_enthalpie = 433107
    benzene.freezing_T = 278.65
    benzene.compute_molar_volume()
    mix.add_component(benzene)

    c1benzene = ou.component('c1 benzene',list_proportion[7])
    c1benzene.density = 875.7
    c1benzene.molar_weight = 0.09215
    c1benzene.boiling_T = 383.75
    c1benzene.partial_P = 3800
    c1benzene.ref_T_Clau = 298.15
    c1benzene.vap_enthalpie = 412480
    c1benzene.solubility = 0.11
    c1benzene.freezing_T = 178.15
    c1benzene.compute_molar_volume()
    mix.add_component(c1benzene)

    c2benzene = ou.component('c2 benzene',list_proportion[8])
    c2benzene.density = rand(870,874.3)
    c2benzene.molar_weight = 0.106
    c2benzene.boiling_T = rand(409.35,413.15)
    c2benzene.partial_P = rand(1070,1270)
    c2benzene.ref_T_Clau = 298.15
    c2benzene.vap_enthalpie = rand(397853,401733)
    c2benzene.solubility = rand(0.1,0.16)
    c2benzene.freezing_T = rand(178.15, 253.18)
    c2benzene.compute_molar_volume()
    mix.add_component(c2benzene)

    c3benzene = ou.component('c3 benzene',list_proportion[9])
    c3benzene.density = rand(859.3,880)
    c3benzene.molar_weight = 0.12
    c3benzene.boiling_T = rand(425.56,442.65)
    c3benzene.partial_P = rand(256,600)
    c3benzene.ref_T_Clau = 298.15
    c3benzene.vap_enthalpie = rand(375485,398785)
    c3benzene.solubility = rand(0.04,0.055)
    c3benzene.freezing_T = rand(173.55, 235)
    c3benzene.compute_molar_volume()
    mix.add_component(c3benzene)

    c45benzene = ou.component('c4-c5 benzene',list_proportion[10])
    c45benzene.density = rand(860,886)
    c45benzene.molar_weight = rand(134,148.24)
    c45benzene.boiling_T = rand(456.4,478.15)
    c45benzene.solubility = 15.4e-3
    c45benzene.partial_P = 0#supposed to be 0...
    c45benzene.solubility = rand(3.37e-3,30.9e-3)
    c45benzene.freezing_T = rand(185.2, 289.5)
    c45benzene.compute_molar_volume()
    mix.add_component(c45benzene)

    #Next use alcanes
    c10 = ou.component('c10',list_proportion[11])
    c10.density = 730
    c10.molar_weight = 0.142
    c10.boiling_T = 447
    c10.partial_P = 170
    c10.ref_T_Clau = 298.15
    c10.vap_enthalpie = 362676
    c10.solubility = 8.7e-5
    c10.freezing_T = 243.4
    c10.compute_molar_volume()
    mix.add_component(c10)

    c1112 = ou.component('c11-c12',list_proportion[12])
    c1112.density = rand(740,749.5)
    c1112.molar_weight = rand(0.156,0.17)
    c1112.boiling_T = rand(468,489)
    c1112.partial_P = rand(18,55)
    c1112.solubility = rand(4.91e-6,3.5e-4)
    c1112.ref_T_Clau = 298.15
    c1112.freezing_T = rand(247.5, 263.5)
    c1112.compute_enthalpy()
    c1112.compute_molar_volume()
    mix.add_component(c1112)

    c1314 = ou.component('c13-c14',list_proportion[13])
    c1314.density = rand(756,762)
    c1314.molar_weight = rand(0.184,0.198)
    c1314.boiling_T = rand(507,528)
    c1314.solubility = rand(5.5e-5,9.6e-5)
    c1314.partial_P = 0#supposed to be 0...
    c1314.freezing_T = rand(268,278)
    c1314.compute_molar_volume()
    mix.add_component(c1314)

    c1516 = ou.component('c15-c16',list_proportion[14])
    c1516.density = rand(769,773)
    c1516.molar_weight = rand(0.21241,0.22643)
    c1516.boiling_T = rand(540.15,543.15)
    c1516.solubility = 3.6e-5
    c1516.partial_P = 0#supposed to be 0...
    c1516.freezing_T = rand(286.5, 291.33)
    c1516.compute_molar_volume()
    mix.add_component(c1516)

    c1718 = ou.component('c17-c18',list_proportion[15])
    c1718.density = rand(777,778)
    c1718.molar_weight = rand(0.240,0.254)
    c1718.boiling_T = rand(576.5,589.5)
    c1718.partial_P = 0#supposed to be 0...
    c1718.solubility = 0 #put to 0
    c1718.freezing_T = rand(295.1, 302)
    c1718.compute_molar_volume()
    mix.add_component(c1718)

    c1920 = ou.component('c19-c20',list_proportion[16])
    c1920.density = rand(785,789)
    c1920.molar_weight = rand(0.268,0.282)
    c1920.boiling_T = rand(603.5,617.5)
    c1920.partial_P = 0#supposed to be 0...
    c1920.solubility = 0 #put to 0
    c1920.freezing_T = rand(305.15, 310.15)
    c1920.compute_molar_volume()
    mix.add_component(c1920)

    c2125 = ou.component('c21-c25',list_proportion[17])
    c2125.density = rand(792,801)
    c2125.molar_weight = rand(0.297,0.353)
    c2125.boiling_T = rand(632.15,675.15)
    c2125.partial_P = 0#supposed to be 0...
    c2125.solubility = 0 #put to 0
    c2125.freezing_T = rand(313.15, 326.15)
    c2125.compute_molar_volume()
    mix.add_component(c2125)

    c25p= ou.component('c25p',list_proportion[18])
    c25p.density = 817
    c25p.molar_weight = 0.563
    c25p.boiling_T = 795.15
    c25p.partial_P = 0#supposed to be 0...
    c25p.solubility = 0 #put to 0
    c25p.freezing_T(354.15)
    c25p.compute_molar_volume()
    mix.add_component(c25p)

    naphta1 = ou.component('Light naphta',list_proportion[19])
    naphta1.density = rand(1000,1180)
    naphta1.molar_weight = rand(0.128,0.142)
    naphta1.boiling_T = rand(491.15,515)
    naphta1.partial_P = rand(4.91,7.2)
    naphta1.solubility = 23e-3
    naphta1.ref_T_Clau = 293.15
    naphta1.vap_enthalpie = 404929
    naphta1.freezing_T = rand(251, 353.35)
    naphta1.compute_molar_volume()
    mix.add_component(naphta1)

    naphta2 = ou.component('Heavy naphta',list_proportion[20])
    naphta2.density = 994
    naphta2.molar_weight = rand(0.156,0.170)
    naphta2.boiling_T = rand(532.15,546)
    naphta2.partial_P = 4.2
    naphta2.solubility = 8e-3
    naphta2.ref_T_Clau = 293.15
    naphta2.freezing_T = 364.65
    naphta2.compute_enthalpy()
    naphta2.compute_molar_volume()
    mix.add_component(naphta2)

    PAH1 = ou.component('PAH 1',list_proportion[21])
    PAH1.density = rand(898.7,1280)
    PAH1.molar_weight = rand(0.152,0.184)
    PAH1.boiling_T = rand(528,614.5)
    PAH1.partial_P = rand(0,0.7)
    PAH1.ref_T_Clau = 293.15
    PAH1.solubility = rand(0.044e-3, 4.45e-3)
    PAH1.freezing_T = rand(342.3, 489)
    PAH1.compute_molar_volume()
    mix.add_component(PAH1)

    PAH2 = ou.component('PAH 2',list_proportion[22])
    PAH2.density = rand(1232,1378)
    PAH2.molar_weight = rand(0.202,0.278)
    PAH2.boiling_T = rand(552,773)
    PAH2.partial_P = rand(0,0.7)
    PAH2.solubility = rand(2e-6, 0.265e-3)
    PAH2.freezing_T = rand(383.9,551)
    PAH2.compute_molar_volume()
    mix.add_component(PAH2)


    phenol = ou.component('Phenols',list_proportion[23])
    phenol.density = rand(1011,1070)
    phenol.molar_weight = rand(0.094,0.122)
    phenol.boiling_T = rand(454.95,491)
    phenol.partial_P = rand(24,46.6)
    phenol.ref_T_Clau = 298.15
    phenol.vap_enthalpie = 614387
    phenol.solubility = rand(26,84)
    phenol.freezing_T = rand(298.85, 318.2)
    phenol.compute_molar_volume()
    mix.add_component(phenol)


    #Go from %massic to %volumic
    array_volume = []

    for comp in mix.list_component:
        array_volume.append(comp.amount/comp.get_density())

    total = sum(array_volume)
    for id in range(0, len(mix.list_component)):
        mix.list_component[id].amount = array_volume[id]/100 * volume_init

    return mix

def rand(min, max):
    return random.random() *(max-min)+min

def createMixImaros(name, volume_init, list_proportion, density,t_density, viscosity, t_viscosity):
    """
    sum of list_proportion (it is the massic proportion)
    !!!!viscosity is in mPas
    """
    mix = ou.mix(name)
    mix.density = density
    mix.density_T = t_density
    if len(viscosity) > 0:
        for visc in viscosity:
             visc = visc/1000
        mix.viscosity = viscosity
        mix.viscosity_T = t_viscosity

    c1c4 = ou.component('c1-c4',list_proportion[0])
    #c1c4.density = 1.6433#gas
    c1c4.density = 543.75#liquid
    c1c4.molar_weight = 0.037
    c1c4.boiling_T = 200.3
    c1c4.partial_P = 1621053
    #c1c4.vap_enthalpie = 510625#TO VERIFY ! at 100K...
    c1c4.solubility = 46.7e-3
    c1c4.freezing_T = 100.5
    c1c4.compute_molar_volume()
    mix.add_component(c1c4)

    c5 = ou.component('c5 sat',list_proportion[1])
    c5.density = 668.45
    c5.molar_weight = 0.0713
    c5.boiling_T = 310.84
    c5.partial_P = 59773
    c5.ref_T_Clau = 293.15
    c5.vap_enthalpie = 366190
    c5.solubility = 0.0915
    c5.freezing_T = 145.2
    c5.compute_molar_volume()
    mix.add_component(c5)

    c6 = ou.component('c6 sat',list_proportion[2])
    c6.density = 700.2
    c6.molar_weight = 0.0855
    c6.boiling_T = 343.34
    c6.partial_P = 17998.5
    c6.ref_T_Clau = 300.65
    c6.vap_enthalpie = 378552.5
    c6.solubility = 25.8e-3
    c6.freezing_T = 192.63
    c6.compute_molar_volume()
    mix.add_component(c6)

    c7 = ou.component('c7 sat',list_proportion[3])
    c7.density = 712.6
    c7.molar_weight = 0.0995
    c7.boiling_T = 367.75
    c7.partial_P = 6133
    c7.ref_T_Clau = 298.15
    c7.vap_enthalpie = 364960
    c7.solubility = 7.24e-3
    c7.freezing_T = 161.15
    c7.compute_molar_volume()
    mix.add_component(c7)

    c8 = ou.component('c8 sat',list_proportion[4])
    c8.density = 698.3
    c8.molar_weight = 0.114
    c8.boiling_T = 395.075
    c8.partial_P = 5300
    c8.ref_T_Clau = 310.15
    c8.solubility = 0.66e-3
    c8.freezing_T = 152
    c8.compute_enthalpy()
    c8.compute_molar_volume()
    mix.add_component(c8)

    c9 = ou.component('c9 sat',list_proportion[5])
    c9.density = 718
    c9.molar_weight = 0.128
    c9.boiling_T = 424
    c9.partial_P = 590
    c9.ref_T_Clau = 298.15
    c9.solubility = 0 #put to 0
    c9.freezing_T = 219.5
    c9.compute_enthalpy()
    c9.compute_molar_volume()
    mix.add_component(c9)

    benzene = ou.component('benzene',list_proportion[6])
    benzene.density = 880
    benzene.molar_weight = 0.07811
    benzene.boiling_T = 353.25
    benzene.partial_P = 12700
    benzene.ref_T_Clau = 298.15
    benzene.solubility = 1.360
    benzene.vap_enthalpie = 433107
    benzene.freezing_T = 278.65
    benzene.compute_molar_volume()
    mix.add_component(benzene)

    c1benzene = ou.component('c1 benzene',list_proportion[7])
    c1benzene.density = 875.7
    c1benzene.molar_weight = 0.09215
    c1benzene.boiling_T = 383.75
    c1benzene.partial_P = 3800
    c1benzene.ref_T_Clau = 298.15
    c1benzene.vap_enthalpie = 412480
    c1benzene.solubility = 0.11
    c1benzene.freezing_T = 178.15
    c1benzene.compute_molar_volume()
    mix.add_component(c1benzene)

    c2benzene = ou.component('c2 benzene',list_proportion[8])
    c2benzene.density = 872.15
    c2benzene.molar_weight = 0.106
    c2benzene.boiling_T = 411.25
    c2benzene.partial_P = 1170
    c2benzene.ref_T_Clau = 298.15
    c2benzene.vap_enthalpie = 399793
    c2benzene.solubility = 0.13
    c2benzene.freezing_T = 215.665
    c2benzene.compute_molar_volume()
    mix.add_component(c2benzene)

    c3benzene = ou.component('c3 benzene',list_proportion[9])
    c3benzene.density = 867.77
    c3benzene.molar_weight = 0.12
    c3benzene.boiling_T = 433.53
    c3benzene.partial_P = 437
    c3benzene.ref_T_Clau = 298.15
    c3benzene.vap_enthalpie = 386274
    c3benzene.solubility = 0.483
    c3benzene.freezing_T = 195.18
    c3benzene.compute_molar_volume()
    mix.add_component(c3benzene)

    c45benzene = ou.component('c4-c5 benzene',list_proportion[10])
    c45benzene.density = 868.67
    c45benzene.molar_weight = 0.138747
    c45benzene.boiling_T = 468.517
    c45benzene.solubility = 15.4e-3
    c45benzene.partial_P = 0#supposed to be 0...
    c45benzene.solubility = 0 #put to 0
    c45benzene.freezing_T = 237.35
    c45benzene.compute_molar_volume()
    mix.add_component(c45benzene)

    #Next use alcanes
    c10 = ou.component('c10',list_proportion[11])
    c10.density = 730
    c10.molar_weight = 0.142
    c10.boiling_T = 447
    c10.partial_P = 170
    c10.ref_T_Clau = 298.15
    c10.vap_enthalpie = 362676
    c10.solubility = 8.7e-5
    c10.freezing_T = 243.4
    c10.compute_molar_volume()
    mix.add_component(c10)

    c1112 = ou.component('c11-c12',list_proportion[12])
    c1112.density = 744.75
    c1112.molar_weight = 0.163
    c1112.boiling_T = 478.5
    c1112.partial_P = 36.5
    c1112.solubility = 1.7e-4
    c1112.ref_T_Clau = 298.15
    c1112.freezing_T = 255.5
    c1112.compute_enthalpy()
    c1112.compute_molar_volume()
    mix.add_component(c1112)

    c1314 = ou.component('c13-c14',list_proportion[13])
    c1314.density = 759
    c1314.molar_weight = 0.191
    c1314.boiling_T = 517.5
    c1314.solubility = 7.55e-5
    c1314.partial_P = 0#supposed to be 0...
    c1314.freezing_T = 273
    c1314.compute_molar_volume()
    mix.add_component(c1314)

    c1516 = ou.component('c15-c16',list_proportion[14])
    c1516.density = 771
    c1516.molar_weight = 0.219
    c1516.boiling_T = 541.65
    c1516.solubility = 3.6e-5
    c1516.partial_P = 0#supposed to be 0...
    c1516.freezing_T = 288.915
    c1516.compute_molar_volume()
    mix.add_component(c1516)

    c1718 = ou.component('c17-c18',list_proportion[15])
    c1718.density = 777.5
    c1718.molar_weight = 0.247
    c1718.boiling_T = 583
    c1718.partial_P = 0#supposed to be 0...
    c1718.solubility = 0 #put to 0
    c1718.freezing_T = 298.55
    c1718.compute_molar_volume()
    mix.add_component(c1718)

    c1920 = ou.component('c19-c20',list_proportion[16])
    c1920.density = 787
    c1920.molar_weight = 0.275
    c1920.boiling_T = 610.5
    c1920.partial_P = 0#supposed to be 0...
    c1920.solubility = 0 #put to 0
    c1920.freezing_T = 307.65
    c1920.compute_molar_volume()
    mix.add_component(c1920)

    c2125 = ou.component('c21-c25',list_proportion[17])
    c2125.density = 793
    c2125.molar_weight = 0.325
    c2125.boiling_T = 653.55
    c2125.partial_P = 0#supposed to be 0...
    c2125.solubility = 0 #put to 0
    c2125.freezing_T = 319.95
    c2125.compute_molar_volume()
    mix.add_component(c2125)

    c25p= ou.component('c25p',list_proportion[18])
    c25p.density = 984
    c25p.molar_weight = 0.563
    c25p.boiling_T = 795.15
    c25p.partial_P = 0#supposed to be 0...
    c25p.solubility = 0 #put to 0
    c25p.freezing_T = 354.15
    c25p.compute_molar_volume()
    mix.add_component(c25p)

    naphta1 = ou.component('Light naphta',list_proportion[19])
    naphta1.density = 1090
    naphta1.molar_weight = 0.135
    naphta1.boiling_T = 503
    naphta1.partial_P = 6
    naphta1.solubility = 23e-3
    naphta1.ref_T_Clau = 293.15
    naphta1.vap_enthalpie = 404929
    naphta1.freezing_T = 302.175
    naphta1.compute_molar_volume()
    mix.add_component(naphta1)

    naphta2 = ou.component('Heavy naphta',list_proportion[20])
    naphta2.density = 994
    naphta2.molar_weight = 0.163
    naphta2.boiling_T = 539.075
    naphta2.partial_P = 4.2
    naphta2.solubility = 8e-3
    naphta2.ref_T_Clau = 293.15
    naphta2.freezing_T = 365.65
    naphta2.compute_enthalpy()
    naphta2.compute_molar_volume()
    mix.add_component(naphta2)

    PAH1 = ou.component('PAH 1',list_proportion[21])
    PAH1.density = 1125
    PAH1.molar_weight = 0.165
    PAH1.boiling_T = 575
    PAH1.partial_P = 0.1
    PAH1.ref_T_Clau = 293.15
    PAH1.solubility = 1.05e-3
    PAH1.freezing_T = 385.4
    PAH1.compute_molar_volume()
    mix.add_component(PAH1)

    PAH2 = ou.component('PAH 2',list_proportion[21])
    PAH2.density = 1277
    PAH2.molar_weight = 0.247
    PAH2.boiling_T = 689.8
    PAH2.solubility = 0.139e-3
    PAH2.partial_P = 0
    PAH2.freezing_T = 488.625
    PAH2.compute_molar_volume()
    mix.add_component(PAH2)

    phenol = ou.component('Phenols',list_proportion[23])
    phenol.density = 1037
    phenol.molar_weight = 0.108
    phenol.boiling_T = 472.3
    phenol.partial_P = 35.3
    phenol.ref_T_Clau = 298.15
    phenol.vap_enthalpie = 614387
    phenol.solubility = 55
    phenol.freezing_T = 311.07
    phenol.compute_molar_volume()
    mix.add_component(phenol)


    #Go from %massic to %volumic
    array_volume = []

    for comp in mix.list_component:
        array_volume.append(comp.amount/comp.get_density())

    total = sum(array_volume)
    for id in range(0, len(mix.list_component)):
        mix.list_component[id].amount = array_volume[id]/total * volume_init

    return mix

def createMixDisti(path, name, density, t_density, viscosity, t_viscosity, tot_amount):
    array_T = []
    array_frac = []
    with open(path+"_lab.csv", newline='') as csvfile:
     csvfile = csv.reader(csvfile, delimiter=',', quotechar='"')
     for row in csvfile:
         if(len(array_T) > 0):
             if(array_T[-1] < int(row[0]) and array_frac[-1] < float(row[1])):
                 array_T.append(int(row[0]))
                 array_frac.append(float(row[1]))
         else:#the first element should be always added
             array_T.append(int(row[0]))
             array_frac.append(float(row[1]))
    #simulation
    with open(path+"_sim.csv", newline='') as csvfile:
     csvfile = csv.reader(csvfile, delimiter=',', quotechar='"')
     for row in csvfile:
         if(array_T[-1] < int(row[0]) and array_frac[-1] < float(row[1])):
             array_T.append(int(row[0]))
             array_frac.append(float(row[1]))

    mix = ou.mix(name)
    mix.density = density
    mix.density_T = t_density
    if len(viscosity) > 0:
        for visc in viscosity:
             visc = visc/1000
        mix.viscosity = viscosity
        mix.viscosity_T = t_viscosity
    mix.generate_component_cut(array_T, array_frac, tot_amount)
    mix.add_oil_properties()
    return mix


def load_ADIOS(path, oilName, oil_amount):
    with open(path, newline='') as csvfile:
     csvfile = csv.reader(csvfile, delimiter=';', quotechar='"')

     for row in csvfile:
         if row[0] ==  oilName:
             return build_ADIOS(row, oil_amount)

     print('oil not found')

def build_ADIOS(row, oil_amount):
    cut_T = []
    cut_fract = []
    for i in range(0,15):
        fract = row[73+i*3]
        temp = row[71+i*3]
        if fract != '' and temp != '':
            temp = temp.replace(',', '.')
            fract = fract.replace(',', '.')
            cut_T.append(float(temp)-273.15)
            cut_fract.append(float(fract)*100)
    mix = ou.mix(row[0])
    mix.generate_component_cut(cut_T, cut_fract, oil_amount)

    density = []
    density_T = []
    for i in range(0,3):
        dens = row[23+i*3]
        temp = row[24+i*3]
        if dens != '' and temp != '':
            density_T.append(float(temp))
            density.append(float(dens))

    if len(density) > 0:
        mix.density = density
        mix.density_T = density_T
    return mix

def compare_oil_properties(oil, T):
    """
    Compute the boiling temperature, molar volume and molar weigth of
    components

    Parameters
    ----------
    oil : Oil mix
    T : Temperature [K]

    """
    for component in oil.list_component:
        #ebulition temperature
        boiling_T = ev.boiling_T_rho(component.density)
        partial_P = ev.vapor_pressure_eb_T(component.boiling_T,T)
        molar_volume = ev.molar_volume_eb_T(component.boiling_T)
        if(component.partial_P is None):
            component.partial_P = 0

        print(component.name+"-----------")
        print("boiling T° [K] real:"+str(int(component.boiling_T))+"  computed:"+str(int(boiling_T))+"          % :"+str(int(100*component.boiling_T/boiling_T)))
        print("partial_P [Pa] real:"+str(int(component.partial_P))+"  computed:"+str(int(partial_P))+"          % :"+str(int(100*component.partial_P/partial_P)))
        print("molar_volume [l/kmol] real:"+str(int(component.molar_volume*1000000))+"  computed:"+str(int(molar_volume*1000000))+"          % :"+str(int(100*component.molar_volume/molar_volume)))
        print("  ")

def saveOil(volume_init):
    """
    Save the IMAROS oil in json
    """
    #IM1  = createMixImaros("to_delete" , volume_init, [0.74,1,0.01,1,1,0.01,1,1,0.01,0.01,1,0.01,2.76,3.01,3.01,3.02,4.46,19.30,61.16,1.36,0.30,0.52,0.35,0], [0.96, 0.95], [278.15, 288.15], [], [])

    IM1  = createMixImaros("IM-1" , volume_init, [0.74,0,0.00,0,0,0.00,0,0,0.00,0.01,0,0.00,2.76,3.01,3.01,3.02,4.46,19.30,61.16,1.36,0.30,0.52,0.35,0], [0.96, 0.95], [278.15, 288.15], [], [])
    IM2  = createMixImaros("IM-2" , volume_init, [0.00,0,0.01,0,0,0.00,0,0,0.00,0.00,0,0.04,0.03,0.03,0.12,0.48,1.42,10.89,86.72,0.00,0.01,0.04,0.21,0], [0.94, 0.93], [278.15, 288.15], [], [])
    IM3  = createMixImaros("IM-3" , volume_init, [0.34,0,0.01,0,0,0.00,0,0,0.00,0.01,0,0.06,3.32,5.12,5.54,5.60,5.37,12.94,57.54,1.34,0.23,1.57,1.01,0], [0.99, 0.98], [278.15, 288.15], [4858, 1293], [278.15, 288.15])
    IM4  = createMixImaros("IM-4" , volume_init, [0.03,0,0.00,0,0,0.00,0,0,0.01,0.02,0,0.58,3.40,5.49,8.14,6.80,6.30,11.43,52.73,1.18,1.05,1.44,1.40,0], [0.95, 0.95], [278.15, 288.15], [2808, 703], [278.15, 288.15])
    IM5  = createMixImaros("IM-5" , volume_init, [0.08,0,0.00,0,0,0.02,0,0,0.00,0.02,0,2.09,3.88,5.32,5.31,4.97,4.04,08.86,63.02,0.17,0.22,0.29,0.54,0], [0.92, 0.91], [278.15, 288.15], [1826, 375], [278.15, 288.15])
    IM6  = createMixImaros("IM-6" , volume_init, [0.00,0,0.00,0,0,0.00,0,0,0.00,0.22,0,1.92,7.52,6.43,6.05,4.60,3.86,05.82,58.65,2.78,0.30,1.31,0.54,0], [0.98, 0.97], [278.15, 288.15], [2244, 892], [278.15, 288.15])
    IM7  = createMixImaros("IM-7" , volume_init, [0.17,0,0.00,0,0,0.00,0,0,0.00,0.00,0,0.32,2.15,3.66,5.73,5.61,5.83,13.22,61.45,0.37,0.20,0.39,0.90,0], [0.95, 0.94], [278.15, 288.15], [4415, 19117], [278.15, 288.15])
    IM8  = createMixImaros("IM-8" , volume_init, [0.53,0,0.00,0,0,0.00,0,0,0.00,0.00,0,0.00,3.62,2.88,2.98,2.76,3.63,11.04,69.38,1.08,0.22,0.54,1.33,0], [0.97, 0.96], [278.15, 288.15], [15585, 3348], [278.15, 288.15])
    IM9  = createMixImaros("IM-9" , volume_init, [1.06,0,0.00,0,0,0.00,0,0,0.00,0.00,0,0.83,5.29,5.41,4.99,4.97,4.55,10.73,61.04,0.58,0.11,0.19,0.25,0], [0.90, 0.90], [278.15, 288.15], [], [])
    IM10 = createMixImaros("IM-10", volume_init, [8.52,0,0.00,0,0,0.00,0,0,0.00,0.01,0,0.53,0.77,1.68,2.39,2.80,3.86,14.84,63.32,0.05,0.16,0.33,0.74,0], [0.95, 0.94], [278.15, 288.15], [12443, 2451], [278.15, 288.15])
    IM11 = createMixImaros("IM-11", volume_init, [0.06,0,0.00,0,0,0.00,0,0,0.00,0.01,0,0.57,0.84,1.84,2.63,2.94,4.08,16.59,69.08,0.06,0.18,0.37,0.75,0], [0.95, 0.94], [278.15, 288.15], [8171, 1964], [278.15, 288.15])
    IM12 = createMixImaros("IM-12", volume_init, [0.00,0,0.00,0,0,0.01,0,0,0.00,0.02,0,1.91,6.24,6.00,5.70,4.72,4.44,07.00,62.04,0.51,0.23,0.33,0.85,0], [0.95, 0.94], [278.15, 288.15], [10679, 3042], [278.15, 288.15])
    IM13 = createMixImaros("IM-13", volume_init, [0.00,0,0.00,0,0,0.01,0,0,0.00,0.01,0,2.65,5.01,4.62,3.07,2.28,2.72,05.99,72.04,0.42,0.14,0.24,0.73,0], [0.96, 0.96], [278.15, 288.15], [24994, 6240], [278.15, 288.15])
    IM14 = createMixImaros("IM-14", volume_init, [0.07,0,0.00,0,0,0.00,0,0,0.00,0.00,0,0.00,1.11,0.92,1.12,1.35,2.36,09.65,82.74,0.02,0.10,0.14,0.42,0], [0.945, 0.937], [278.15, 288.15], [21007, 17121], [278.15, 288.15])
    IM15 = createMixImaros("IM-15", volume_init, [0.03,0,0.00,0,0,0.01,0,0,0.00,0.02,0,0.00,4.01,4.36,4.47,3.29,3.33,05.94,71.66,0.64,0.40,0.46,1.38,0], [0.958, 0.951], [278.15, 288.15], [19406, 4305], [278.15, 288.15])

    list_oil = []
    list_oil.append(IM1)
    list_oil.append(IM2)
    list_oil.append(IM3)
    list_oil.append(IM4)
    list_oil.append(IM5)
    list_oil.append(IM6)
    list_oil.append(IM7)
    list_oil.append(IM8)
    list_oil.append(IM9)
    list_oil.append(IM10)
    list_oil.append(IM11)
    list_oil.append(IM12)
    list_oil.append(IM13)

    for oil in list_oil:
        iojson.save_mix(oil,"oil_oscar/"+oil.name+".json")



volume_init = 0.02
saveOil(volume_init)
depth = 0.9
water_volume = 7
surface = 7.78
wind_speed = 5
current_speed = 0.4
wave_height =  0.75
temperature = 5+273.15
stability_class = 'C'



dt = 60
sim_length = 3600*24*7

surface_array = np.ones(int(sim_length / dt))

time_to_full_area = 0*3600/dt #Instanly on all the surface
for i in range(0, len(surface_array)):
    if i < time_to_full_area:
        surface_array[i] = surface * ((i+0.1)/time_to_full_area)
    else :
        surface_array[i] = surface


#pseudocomponent
IM5  = createMixImaros("IM-5" , volume_init, [0.08,0,0.00,0,0,0.02,0,0,0.00,0.02,0,2.09,3.88,5.32,5.31,4.97,4.04,08.86,63.02,0.17,0.22,0.29,0.54,0], [0.92, 0.91], [278.15, 288.15], [3051, 507], [278.15, 288.15])
IM14 = createMixImaros("IM-14", volume_init, [0.07,0,0.00,0,0,0.00,0,0,0.00,0.00,0,0.00,1.11,0.92,1.12,1.35,2.36,09.65,82.74,0.02,0.10,0.14,0.42,0], [0.945, 0.937], [278.15, 288.15], [71747, 17121], [278.15, 288.15])
IM15 = createMixImaros("IM-15", volume_init, [0.03,0,0.00,0,0,0.01,0,0,0.00,0.02,0,0.00,4.01,4.36,4.47,3.29,3.33,05.94,71.66,0.64,0.40,0.46,1.38,0], [0.958, 0.951], [278.15, 288.15], [19406, 4305], [278.15, 288.15])

#DISTILLATION
#IM5 = createMixDisti('oil_disti/IM5', 'IM-5',[0.92, 0.91], [278.15, 288.15], [3051, 507], [278.15, 288.15], volume_init)
#IM14 = createMixDisti('oil_disti/IM14', 'IM-14',[0.945, 0.937], [278.15, 288.15], [71747, 17121], [278.15, 288.15], volume_init)
#IM15 = createMixDisti('oil_disti/IM15', 'IM-15',[0.958, 0.951], [278.15, 288.15], [19406, 4305], [278.15, 288.15], volume_init)




#5°C
IM5.K_em = 17.56
IM14.K_em = 199.65
IM15.K_em = 56.41
#labo
IM5.max_water = 0.71
IM14.max_water = 0.0
IM15.max_water = 0.50
#polludrome average après 20h
#IM5.max_water = 0.8065
#IM14.max_water = 0.171
#IM15.max_water = 0.631


#15°C
#IM5.K_em = 15.99
#IM14.K_em = 18.71
#IM15.K_em = 43.59
#labo
#IM5.max_water = 0.81
#IM14.max_water = 0.5
#IM15.max_water = 0.70
#polludrome average après 20h
#IM5.max_water = 0.8558
#IM14.max_water = 0.566667
#IM15.max_water = 0.7119

mat = wu.compute_weathering(copy.deepcopy(IM5), temperature, wind_speed, sim_length, dt, water_volume,
                None, surface_array,apply_evaporation = 2, apply_dissolution = 1, apply_emulsion = 1, apply_volatilization = 0, wave_height = wave_height, stability_class=stability_class, current_speed = current_speed)
wu.create_csv(copy.deepcopy(IM5), mat, "IM-5")
wu.plot_matrix_mix(IM5, mat)

mat = wu.compute_weathering(copy.deepcopy(IM14), temperature, wind_speed, sim_length, dt, water_volume,
                None, surface_array,apply_evaporation = 2, apply_dissolution = 1, apply_emulsion = 1, apply_volatilization = 0, wave_height = wave_height, stability_class=stability_class, current_speed = current_speed)
wu.create_csv(copy.deepcopy(IM14), mat, "IM-14")
wu.plot_matrix_mix(IM14, mat)

mat = wu.compute_weathering(copy.deepcopy(IM15), temperature, wind_speed, sim_length, dt, water_volume,
                None, surface_array,apply_evaporation = 2, apply_dissolution = 1, apply_emulsion = 1, apply_volatilization = 0, wave_height = wave_height, stability_class=stability_class, current_speed = current_speed)
wu.create_csv(copy.deepcopy(IM15), mat, "IM-15")
wu.plot_matrix_mix(IM15, mat)
