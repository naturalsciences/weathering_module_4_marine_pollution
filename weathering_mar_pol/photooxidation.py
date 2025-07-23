# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 14:27:28 2021

@author: Ludovic Lepers
"""
import math

def simple_half_live(q, k, dt):
    """
    Return the remaining quantity after photooxidation occuring during a time interval

    Parameters
    ----------
    q : Amount
    k : Half life constant [1/s]
    dt : Time interval [s]


    """
    return q*(1-math.e**(-k*dt))
