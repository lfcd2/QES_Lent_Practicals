from P4main import ocean_model_p4
import numpy as np
from OceanTools.tools import plot
import matplotlib.pyplot as plt
# global variables
V_ocean = 1.34e18  # volume of the ocean in m3
SA_ocean = 358e12  # surface area of the ocean in m2
fSA_hilat = 0.15  # fraction of ocean surface area in 'high latitude' box

# variables used to calculate Q
Q_alpha = 1e-4
Q_beta = 7e-4
Q_k = 8.3e17

# salinity balance - the total amount of salt added or removed to the surface boxes
Fw = 0.1  # low latitude evaporation - precipitation in units of m yr-1
Sref = 35  # reference salinity in units of g kg-1

E = Fw * SA_ocean * (1 - fSA_hilat) * Sref  # amount of salt removed from the low lat box, g kg-1 yr-1, ~ kg m-3 yr-1


def initialise_dicts():
    """Calculates the initial values for each box, and saves them in a dictionary

    Returns
    -------
    init_lolat : dict
        initial values for lolat
    init_hilat : dict
        initial values for hilat
    init_deep : dict
        initial values for deep
    init_atmos : dict
        initial values for atmos
    """

    # set up boxes with all the required variables
    # High Latitude Box
    init_hilat = {
        'name': 'hilat',
        'depth': 200,  # box depth, m
        'SA': SA_ocean * fSA_hilat,  # box surface area, m2
        'T': 3.897678,  # initial water temperature, Celcius
        'S': 34.37786,  # initial salinity
        'T_atmos': 0.,  # air temperature, Celcius
        'tau_M': 100.,  # timescale of surface-deep mixing, yr
        'tau_T': 2.,  # timescale of temperature exchange with atmosphere, yr
        'E': -E,  # salt added due to evaporation - precipitation, kg m-3 yr-1
        # Add new variables here
        'DIC': 2.300182,
        'TA': 2.3134328,
        'tau_CO2': 2.
    }
    init_hilat['V'] = init_hilat['SA'] * init_hilat['depth']  # box volume, m3

    # Low Latitude Box
    init_lolat = {
        'name': 'lolat',
        'depth': 100,  # box depth, m
        'SA': SA_ocean * (1 - fSA_hilat),  # box surface area, m2
        'T': 23.60040,  # initial water temperature, Celcius
        'S': 35.37898,  # initial salinity
        'T_atmos': 25.,  # air temperature, Celcius
        'tau_M': 250.,  # timescale of surface-deep mixing, yr
        'tau_T': 2.,  # timescale of temperature exchange with atmosphere, yr
        'E': E,  # salinity balance, PSU m3 yr-1
        'DIC': 2.23882,
        'TA': 2.313432,
        'tau_CO2': 2.
    }
    init_lolat['V'] = init_lolat['SA'] * init_lolat['depth']  # box volume, m3

    # Deep Ocean Box
    init_deep = {
        'name': 'deep',
        'V': V_ocean - init_lolat['V'] - init_hilat['V'],  # box volume, m3
        'T': 5.483637,  # initial water temperature, Celcius
        'S': 34.47283,  # initial salinity
        # Add new variables here
        'DIC': 2.295687,
        'TA': 2.313432,
    }

    # Atmosphere box
    init_atmos = {
        'name': 'atmos',
        'mass': 5e21,  # in grams
        'moles_air': 1.736e20,
        'moles_CO2': 2.212944e+17,
        'GtC_emissions': 0
    }
    init_atmos['pCO2'] = 1e6 * init_atmos['moles_CO2'] / init_atmos['moles_air']

    # returns the dictionaries
    return init_lolat, init_hilat, init_deep, init_atmos


def question_4():

    dicts = list(initialise_dicts())

    # create a new time axis for the model containing 3000 years with a 0.5 year time step
    tmax = 3000  # how many years to simulate (yr)
    dt = 0.5  # the time step of the simulation (yr)
    time = np.arange(0, tmax + dt, dt)  # the time axis for the model

    # create an array containing GtC_emissions that contains zeros except between years 800-1000
    # where 8 GtC are emitted each year.
    emit_atmos = dicts[-1].copy()  # create a copy of the original atmosphere input dictionary
    emit_atmos['GtC_emissions'] = np.zeros(time.shape)  # creat an array to hold the emission scenario
    emit_atmos['GtC_emissions'][(time > 500) & (time <= 700)] = 8.0  # set e to 8 GtC per year between 500-700
    dicts[-1] = emit_atmos

    # run the model using this emission scenario, and create the required plot
    time_array, finished_dicts = ocean_model_p4(dicts, tmax, dt)
    final_lolat, final_hilat, final_deep, final_atmos = finished_dicts
    fig, axs = plot.boxes(time_array, ['DIC', 'pCO2', 'GtC_emissions'],
                          final_lolat, final_hilat, final_deep, final_atmos)

    plt.show()


if __name__ == "__main__":
    question_4()
