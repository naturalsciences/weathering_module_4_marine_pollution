# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 13:20:25 2021

@author: Ludovic Lepers
"""

import math
import numpy as np
import evaporation as ev
import weathering_utils as wu
import viscosity as vis

def interp(val, array_value, array_ref):
    """
    Interpolate the array_value from the array_ref with val. The array_ref
    must be in an increasing order!
    """
    if val <= array_ref[0]:
        return array_value[0]
    elif val > array_ref[len(array_ref)-1]:
        return array_value[len(array_ref)-1]
    else:
        i = 1
        while val > array_ref[i]:
            i += 1
        delta = array_ref[i] - array_ref[i-1]
        return ((array_ref[i]-val)*array_value[i-1]+array_value[i]*(val-array_ref[i-1]))/delta

class mix:

    def __init__(self, name):
        self.name = name
        self.list_component = []
        self.density = [] #density values depending on the temperature
        self.density_T = [] #temperature of the density values given above
        self.viscosity = [] #density values depending on the temperature
        self.viscosity_T = [] #temperature of the density values given above
        self.K_em = None
        self.max_water = 0.8 # max amount of water in the emulsion


    def generate_component_cut(self, temp, fraction, tot_amount):
        """
        Generate the components for oils, two vectors (one with temperatures
        and one with cumulative amount evaporated). The total amount (m³) must
        be given and the temperature at which the fraction does not evaporated
        anymore too (in °C). It compute some basic properties for the components


        Parameters
        ----------
        temp : Temperature vector
        fraction : Fraction vector
        tot_amount : Amount of the mix [m³]

        Raises
        ------
        Exception
            Raise an exception if the two vectors does not have the same size

        """


        if len(temp) != len(fraction):
            raise Exception("The two vectors must have the same size")

        #Temperature is given in °C
        prev = 0
        remaining = 100
        for i in range(len(temp)):
            if temp[i] > wu.MAX_EVAPORATIVE_TEMP-273.15:
                break
            ratio = fraction[i] - prev
            prev = fraction[i]
            compound = component(temp[i], tot_amount * ratio/100)
            compound.boiling_T = temp[i] + 273.15
            self.list_component.append(compound)
            remaining -= ratio

        #add non volatils compounds if not = 100%
        if remaining > 0:
            compound = component("Heavy (1000)", tot_amount * remaining/100)
            compound.boiling_T = 1000
            self.list_component.append(compound)

        self.add_oil_properties()

    def add_Fingas(self, c1, c2 = None):
        """
        Add the fingas constant c1 and c2 to the mix

        """
        self.fingas1 = c1
        self.fingas2 = c2


    def add_component(self, compound):
        """
        Add a component to the mix

        Parameters
        ----------
        compound : Component to add

        """
        self.list_component.append(compound)

    def is_pure(self):
        """
        Return True if the mix is only composed of 0 or one component, else
        return False.

        """
        return len(self.list_component) <= 1

    def get_comp(self, index):
        """
        Return the component at the position 'index' in the mix

        """
        return self.list_component[index]

    def get_prop(self, T):
        """
        Return the properties of the mix. If the mix is empty, return an
        exception. If it contains only one component, it returns it. If it
        contains more, it returns a combinaison of the properties

        Parameters
        ----------
        T: temperature [K]

        """
        if len(self.list_component) <= 0:
            raise Exception("This mixture has no component!")
        elif len(self.list_component) == 1:
            return self.list_component[0]
        else:   #combine the components
            tot_amount = 0
            molar_sum = 0
            mass_sum = 0
            for comp in self.list_component:
                tot_amount += comp.amount
                molar_sum += comp.amount / comp.molar_volume
                if comp.density is not None:
                    mass_sum +=  comp.amount*comp.density
                else:
                    mass_sum = None

            partial_P = 0
            self.ref_T_Clau = self.list_component[0].ref_T_Clau
            for comp in self.list_component:
                if comp.partial_P is not None:
                    partial_P += comp.get_partial_P(T) * comp.amount / comp.molar_volume
                else:
                    partial_P += 0

            if partial_P is not None and tot_amount != 0:
                partial_P = partial_P / molar_sum
            else:
                partial_P = None


            #1: computing the density from the components
            density = 0
            for comp in self.list_component:
                if comp.get_density(T) is not None and comp.amount > 0:
                    density += comp.amount * comp.get_density(T) / tot_amount

            #2: if not from components, take the defaut
            if density == 0:
                if len(self.density) > 0:
                    density = self.get_density(T)
                elif mass_sum is not None and tot_amount != 0:
                    density = mass_sum/tot_amount
                else:
                    density = None



            bulk = component('bulk', tot_amount)

            bulk.density = density
            if tot_amount != 0:
                bulk.molar_weight = (tot_amount*density)/ molar_sum

                bulk.molar_volume = tot_amount /molar_sum
            else:
                bulk.molar_weight = None
                bulk.molar_volume = None

            bulk.partial_P = partial_P

            return bulk



    def get_molar_fract(self, tg_comp):
        """
        Get the molar fraction of the component tg_comp

        """
        if len(self.list_component) <= 0:
            raise Exception("This oil has no component!")

        molar_sum = 0

        for comp in self.list_component:
            molar_sum += comp.amount / comp.molar_volume

        fract = (tg_comp.amount / tg_comp.molar_volume) /molar_sum
        return fract


    def add_amount(self, add_amount):
        """
        This function add/remove for each component at the same time and return
        an array with the quantities added to each component
        """

        tot_amount = 0
        array_tr = np.zeros((len(self.list_component)))
        for comp in self.list_component:
            tot_amount += comp.amount

        for i in range(0, len(self.list_component)):
            amount = add_amount * (self.list_component[i].amount/tot_amount)
            self.list_component[i].amount += amount
            array_tr[i] = amount
        return array_tr


    def get_density(self, T):
        """
        Return the density[kg/m³] interpolated at the value T [K]
        """
        return interp(T, self.density, self.density_T)

    def get_viscosity(self, T):
        """
        Return the viscosity[Pa s] interpolated at the value T [K]
        """
        return interp(T, self.viscosity, self.viscosity_T)


    def get_array_amount(self):
        """
        Return an array with the amount of each component [m³]
        """
        array_tr = np.zeros((len(self.list_component)))

        for i in range(0,len(self.list_component)):
            array_tr[i] = self.get_comp(i).amount
        return array_tr

    def add_oil_properties(self):
        """
        Compute the boiling temperature, molar volume and molar weigth of
        components

        """
        for component in self.list_component:
            if component.boiling_T is None and component.density is not None:
                component.boiling_T = ev.boiling_T_rho(component.density)
            elif component.boiling_T is not None and component.density is None:
                component.density = ev.rho_from_boiling_T(component.boiling_T)

            if component.molar_volume is None:
                component.molar_volume = ev.molar_volume_eb_T(component.boiling_T)

            if component.density is not None and  component.molar_weight is None:
               component.molar_weight = component.density * component.molar_volume

    def get_mix_viscosity(self, temperature):
        """
        Return the viscosity of the mix using viscosity and freezing point of each
        component
        """
        freezing_p_distrib = 1 #normalized
        if len(self.list_component) <= 0:
            raise Exception("This mixture has no component!")
        elif len(self.list_component) == 1:
            molar_sum = self.list_component[0].amount / self.list_component[0].molar_volume
        else:   #combine the components
            molar_sum = 0
            for comp in self.list_component:
                molar_sum += comp.amount / comp.molar_volume
        visc_sum = 0
        for comp in self.list_component:
            x = comp.amount / comp.molar_volume / molar_sum
            if comp.freezing_T is not None:
                if comp.freezing_T <= temperature:
                    freezing_p_distrib -= x

            visc_sum = x * math.log(comp.get_viscosity(temperature))
        viscosity_ideal = mat.pow(math.e,visc_sum)

        return 1/(1/viscosity_ideal * math.pow(freezing_p_distrib,1-freezing_p_distrib))



def get_emulsion_density(mix, T, array_in_emulsion, water_density = 1020):
    """
    Return the density by taking into account each component, T is the temperature,
    array_in_emulsion is an array with the amount of each component in emulsion
    """
    water_volume = 0
    oil_volume = 0
    oil_mass = 0
    for i in range(0, len(mix.list_component)):
        component_density = mix.list_component[i].get_density(T)
        if component_density is None :
            return None
        oil_volume += (mix.list_component[i].amount+array_in_emulsion[i])
        oil_mass += (mix.list_component[i].amount+array_in_emulsion[i]) * component_density
        water_volume += array_in_emulsion[i] *(mix.max_water/(1-self.max_water))

    return (oil_mass+water_volume*1020)/(oil_volume+water_volume)


class component:
    ref_T_Clau = None #[K] : ref temp for the vap_enthalpie and partial_P
    molar_weight = None #[kg/mol]
    density = None #[kg/m³]
    boiling_T = None #[K]
    partial_P = None #[Pa]
    molar_volume = None #[m³/mol]
    vap_enthalpie = None #[J/kg]
    h_l_phot = None #[1/s]
    h_l_biod = None #[1/s]
    solubility = None #[kg/m³]
    freezing_T = None #[K]

    def __init__(self, name, amount = 0):
        self.name = name
        self.amount = amount #[m³]

    def get_density(self, T=0):
        """
        Return the density [kg/m³] of the component. The temperature T [K]can be given
        for taking into account the density change with the temperature

        """
        #TODO : change in density with the temperature
        return self.density

    def get_partial_P(self, T):
        """
        Return the partial pressure [Pa] of the component at the temperature T [K]

        """
        if self.vap_enthalpie is not None and self.ref_T_Clau is not None:
            a = (1/self.ref_T_Clau)-(1/T)
            return self.partial_P * math.exp((self.vap_enthalpie *
                                              self.molar_weight/8.314) * a)
        else:
            return self.partial_P

    def compute_enthalpy(self, atm_p = 101325):
        """
        Return an estimation of the vap_enthalpie from a vapor pressure at a
        certain temperature an the ebulition temperature. It suppose the
        vap_enthalpie independant of the temperature

        """
        if self.molar_weight is None:
            raise Exception("The molar_weight must be defined for this")
        if self.boiling_T is None:
            raise Exception("The boiling_T must be defined for this")
        if self.ref_T_Clau is None:
            raise Exception("The ref_T_Clau must be defined for this")

        num = - 8.314 *math.log(self.partial_P/atm_p)
        den = self.molar_weight* (1/self.boiling_T-1/self.ref_T_Clau)
        self.vap_enthalpie =  num/den


    def compute_molar_volume(self):
        """
        Compute the molar volume from the molar weight and density

        """
        self.molar_volume = self.molar_weight/self.density
        #print(self.name, ":", int(self.molar_volume*1000000)/1000," [l/mol]")

    def get_viscosity(self, T):
        """
        Compute the viscosity of the component [Pa s]
        """
        if self.molar_weight is None:
            raise Exception("The molar_weight must be defined for computing viscosity")
        if self.boiling_T is None:
            raise Exception("The boiling_T must be defined for computing viscosity")
        if self.density is None:
            raise Exception("The density must be defined for computing viscosity")
        return vis.viscosity(T,self.boiling_T, self.molar_weight, self.density)

def to_half_life(days):
    """
    Return the constant [1/s] from the half life length [day]

    """
    s= days * 3600*24
    return -math.log(1/2)/s
