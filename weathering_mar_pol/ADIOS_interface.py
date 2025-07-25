from .oil_utils import component
from .oil_utils import mix
import json
import os

def distillation_from_ADIOS(path: str, amount : float):
    """
    Load distillation data from an ADIOS JSON file.

    path is the path to the ADIOS JSON file.
    amount is the volume of oil (m³)
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"No file found at {path}")
    
    with open(path, 'r') as json_file:
        try:
            data = json.load(json_file)
        except(json.JSONDecodeError):
            raise json.JSONDecodeError(f"File {path} is not a valid json file")

    try:
        mix_oil = mix(data["metadata"]["name"])
        

        load_densities(mix_oil, data)
        load_viscosities(mix_oil, data)

        cuts = data["sub_samples"][0]["distillation_data"]["cuts"]

        fracts = []
        temps = []

        for cut in cuts:
            fracts.append(cut["fraction"]["value"])
            temps.append(cut["vapor_temp"]["value"])

        sort_by_temperature(fracts, temps)

        mix_oil.generate_component_cut(temps, fracts, amount)
        mix_oil.add_oil_properties()

    except(KeyError):
        raise KeyError(f"File {path} is not a valid ADIOS file for distillation")
    

    return mix_oil


def OSCAR_from_ADIOS(path : str, properties : dict, amount : float):
    """
    Load oscar characterization data from an ADIOS JSON file.

    path is the path to the ADIOS JSON file.
    properties is a dictionary of the properties of the components
    amount is the volume of oil (m³)
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"No file found at {path}")
    
    with open(path, 'r') as json_file:
        try:
            data = json.load(json_file)
        except(json.JSONDecodeError):
            raise json.JSONDecodeError(f"File {path} is not a valid json file")

    try:
        mix_oil = mix(data["metadata"]["name"])
        

        load_densities(mix_oil, data)
        load_viscosities(mix_oil, data)

        compounds = data["sub_samples"][0]["compounds"]

        for compound in compounds:
            if "OSCAR Characterization" in compound["groups"]:
                name = compound["name"]
                
                if not "value" in compound["measurement"] or not name in properties:
                    continue
                os_comp = component(name, compound["measurement"]["value"] / 100 * amount)
                properties_com = properties[name]
                os_comp.density = properties_com["density"]
                os_comp.molar_weight = properties_com["molar_mass"]
                os_comp.boiling_T = properties_com["boiling_temp"]
                if "solubility" in properties_com:
                    os_comp.solubility = properties_com["solubility"] * os_comp.molar_weight / 1000
                if "partial_pressure" in properties_com:
                    os_comp.partial_P = properties_com["partial_pressure"]
                os_comp.freezing_T = properties_com["freezing_temp"]
                if "vap_enthalpie" in properties_com and "ref_temp" in properties_com:
                    os_comp.ref_T_Clau = properties_com["ref_temp"]
                    os_comp.vap_enthalpie = properties_com["vap_enthalpie"]
                elif "ref_temp" in properties_com:
                    os_comp.ref_T_Clau = properties_com["ref_temp"]
                    os_comp.compute_enthalpy()
                os_comp.compute_molar_volume()

                mix_oil.add_component(os_comp)

    except(KeyError):
        raise KeyError(f"File {path} is not a valid ADIOS file for OSCAR")
    
    return mix_oil


def sort_by_temperature(to_be_sorted : list, temperatures : list):
    combined = sorted(zip(temperatures, to_be_sorted))
    for i, (temp, dens) in enumerate(combined):
        temperatures[i] = temp
        to_be_sorted[i] = dens


def load_densities(mix_oil : mix, data : dict):
    """
    Load the density in the mix from the ADIOS data
    """
    densities = data["sub_samples"][0]["physical_properties"]["densities"]

    for density in densities:
        mix_oil.density.append(density["density"]["value"])
        mix_oil.density_T.append(density["ref_temp"]["value"])
        if density["density"]["unit"] == "g/cm\u00b3":
            mix_oil.density[-1] *= 1000
        else:
            raise Exception(f"Unit {density["density"]["unit"]} not supported")

    sort_by_temperature(mix_oil.density, mix_oil.density_T)

def load_viscosities(mix_oil : mix, data : dict):
    """
    Load the viscosity in the mix from the ADIOS data
    """
    viscosities = data["sub_samples"][0]["physical_properties"]["dynamic_viscosities"]

    for viscosity in viscosities:
        mix_oil.viscosity.append(viscosity["viscosity"]["value"])
        mix_oil.viscosity_T.append(viscosity["ref_temp"]["value"])

        if viscosity["viscosity"]["unit"] == "mPa s":
            mix_oil.viscosity[-1] *= 1000
        else:
            raise Exception(f"Unit {viscosity["viscosity"]["unit"]} not supported")

    sort_by_temperature(mix_oil.viscosity, mix_oil.viscosity_T)