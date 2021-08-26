# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 15:58:16 2021

@author: Ludovic Lepers

Do not hesitate to look into the others python files for more information on the variable dimension
"""

import copy
import weathering_utils as wu  #contains fonctions
import oil_utils as ou  #contains the objects needed

"""
Defining the initials conditions
"""
#amount of each mix
amount_init = 1000 #[m³]
#amount of water in which the component will be dissolved
water_volume = amount_init*1000 #[m³]
#length of a timestep (can be between 0 and 1 or bigger)
dt = 30    #[s]
#length of the simulation
sim_length = 3600 * 24 * 1 #[s]
#speed of the wind at 10 meters
wind_speed = 5 #[m/s]
temperature = 273.15 +20 #[K]
wave_heigth = 1 #[m]
slick_thickness = 2.3e-2 #[m], totally arbitrary here


"""
Creating a mix for an oil and weathering it
"""

brent_blend_cut_T = [40,80,100,120,150,160,180,200,250,300,350,400,500,600,700]
brent_blend_fract = [3.0,4.0,5.0,19.0,22.0,25.0,29.0,32.0,42.0,52.0,62.0,70.0,85.0,95.0,99.0]
brent_blend = ou.mix('BRENT BLEND') #creating a mix for applying weathering on it
#use the two vectors defined earlier to generate pseudo compoenents. The first
#one is about temperature and the second is the cumulative fraction distillated
brent_blend.generate_component_cut(brent_blend_cut_T, brent_blend_fract, amount_init)
#physical caracteristic of the oil
brent_blend.density.append(835)
brent_blend.density_T.append(280)#temperature of the density mesurement
brent_blend.K_em = 20 #for emulsion
brent_blend.visco = 4.5
#fingas constants of the oil, only if used for evaporation
brent_blend.add_Fingas(3.39, 0.048)
#compute some usefull data (molar masse, etc...)
wu.add_oil_properties(brent_blend, temperature)

#!!! The mix will be modified at the end, to avoid that use copy.deepcopy(mix)
mix = copy.deepcopy(brent_blend)

#create a matrix with the state of the mix at each timestep
mat = wu.all_process(mix, temperature, wind_speed, sim_length, dt, water_volume,
                     slick_thickness)

#draw the graph
wu.plot_matrix_mix(mix, mat)


"""
Creating a mix of xyleme and toluene and weathering it

"""

mix2 = ou.mix('Xylene and Toluene')
xyl = ou.component('Xyleme',amount_init/2)#HNS MS 20°C
xyl.density = 870
xyl.molar_volume = 1.220e-4
xyl.molar_weight = 0.10616
xyl.boiling_T = 140.2+273.15
xyl.partial_P = 1070
xyl.ref_T_Clau = 25+273
xyl.vap_enthalpie = 401733
xyl.solubility = 100e-3
mix2.add_component(xyl)


tol = ou.component('Toluene',amount_init/2)#HNS MS 20°C
tol.density = 868.3
tol.molar_volume = 1.061e-4
tol.molar_weight = 0.09215
tol.boiling_T = 110.58+273.15
tol.partial_P = 3800
tol.ref_T_Clau = 25+273
tol.vap_enthalpie = 412480
tol.solubility = 110e-3
#tranform a number of day into a half life
tol.h_l_biod = wu.to_half_life(30)

tol.mu_max = 1.5e-4
tol.ks = 1.96e-3
tol.Y_oil = 1.22


mix2.add_component(tol)

#make the area slowy increase(needs to be in a array the length of the amount of timestep)
slick_thickness = 5e-3 #[m], totally arbitrary here
area_const=amount_init/slick_thickness#[m²]
area = []
for i in range(0, int(sim_length/dt)):
    area.append((i*dt/sim_length)*area_const)

#create a matrix with the state of the mix at each timestep
#the 2 means ALOHA model and the 0 means no emulsification (cf doc of the fonction)
mat = wu.all_process(mix2, temperature, wind_speed, sim_length, dt, water_volume,
                     slick_thickness, area, 'A', 2, 0, wave_heigth)

#draw the graph
wu.plot_matrix_mix(mix2, mat)



"""
Using emulsion with the brent blend

"""
# deepcopy not needed because brent_blend is not used later
mat = wu.all_process(brent_blend, temperature, wind_speed, sim_length, dt, water_volume,
                     slick_thickness, area, 'A', 3, 1, wave_heigth)

#draw the graph
wu.plot_matrix_mix(brent_blend, mat)
