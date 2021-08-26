# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 14:27:04 2021

@author: Ludovic Lepers
"""

import math
import numpy as np

def simple_half_live(q, k, dt):
    """
    Return the photooxidation occuring during a time interval

    Parameters
    ----------
    q : Amount
    k : Half life constant [1/s]
    dt : Time interval [s]


    """
    return q*(1-math.e**(-k*dt))


def biodegradation_gouttes(mix, sim_length, dt, initial_droplet_size, bact_conc, oil_conc):

    srfc_droplet = 4*math.pi*(initial_droplet_size/2)**2
    initial_droplet_vol = 4/3*math.pi*(initial_droplet_size/2)**3
    amount_droplet = mix.get_prop().amount/initial_droplet_vol
    tot_srcf = amount_droplet *  srfc_droplet
    
    time_step_amount = int(sim_length / dt)
    matrix = np.zeros((int(time_step_amount),2))
    matrix[0,1] = mix.get_prop().amount
    for i in range(10,time_step_amount):
        matrix[int(i),0] = i * dt
        comp = mix.get_prop()

        k = (2/comp.density*comp.mu_max/comp.Y_oil*oil_conc/(comp.ks+oil_conc)
             * bact_conc/tot_srcf) 
        x = (1-(1+(k/initial_droplet_size)*i*dt)**1/3)/(i*dt)
        
        # dD=(-2/comp.density*comp.mu_max/comp.Y_oil*oil_conc/(comp.ks+oil_conc)
        #     * bact_conc*dt)
        #r_old = (comp.amount/amount_droplet*(3/4*math.pi))**(1/3)
        #flux = 4/3*math.pi*(r_old**3-(r_old-dD/2)**3)
        # delta = matrix[i,1-1] * x
        # mix.add_amount(-flux) 
        #print(dD)
        matrix[i,1] = matrix[0,1]-x*matrix[0,1]
    return matrix
    
    
    
    
