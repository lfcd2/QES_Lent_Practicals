import os
import sys
import cbsyst
import numpy as np
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
E = Fw * SA_ocean * (1 - fSA_hilat) * Sref  # amount of salt removed from the low latitude box,  g kg-1 yr-1, ~ kg m-3 yr-1

init_hilat = {
    'name': 'hilat',
    'depth': 200,  # box depth, m
    'SA': SA_ocean * fSA_hilat,  # box surface area, m2
    'T': 15.,  # initial water temperature, Celcius
    'S': 34.,  # initial salinity
    'tau_M': 100.,  # timescale of surface-deep mixing, yr
    'T_atmos': 0.,  # air temperature, Celcius
    'tau_T': 2.,  # timescale of temperature exchange with atmosphere, yr
    'E': -E  # salt added due to evaporation - precipitation, kg m-3 yr-1
}
init_hilat['V'] = init_hilat['SA'] *  init_hilat['depth']  # box volume, m3

init_lolat = {
    'name': 'lolat',
    'depth': 100,  # box depth, m
    'SA': SA_ocean * (1- fSA_hilat),  # box surface area, m2
    'T': 15.,  # initial water temperature, Celcius
    'S': 35.,  # initial salinity
    'tau_M': 250.,  # timescale of surface-deep mixing, yr
    'T_atmos': 25.,  # air temperature, Celcius
    'tau_T': 2.,  # timescale of temperature exchange with atmosphere, yr
    'E': -E  # salt added due to evaporation - precipitation, kg m-3 yr-1
}
init_lolat['V'] = init_lolat['SA'] *  init_lolat['depth']  # box volume, m3

init_deep = {
    'name': 'deep',
    'SA': SA_ocean,  # box surface area, m2
    'T': 5.,  # initial water temperature, Celcius
    'S': 34.5  # initial salinity
}
init_hilat['V'] = init_hilat['SA'] *  init_hilat['depth']  # box volume, m3


def ocean_model(lolat, hilat, deep, tmax, dt):
    """Run the ocean model for a given time period and return the results for each box.

    Parameters
    ----------
    lolat, hilat, deep : dict
        dictionaries containing the box properties
    tmax : int or float
        The maximum time to run the model for (yr)
    dt : float
        The time step of the model (yr)

    Returns
    -------
    tuple of (time, lolat, hilat, deep)
    """

    # create the time scale for the model
    time = np.arange(0, tmax + dt, dt)

    # identify which variables will change with time
    model_vars = ['T', 'S']

    # create copies of the input dictionaries so we don't modify the originals
    lolat = lolat.copy()
    hilat = hilat.copy()
    deep = deep.copy()

    # turn all time-evolving variables into arrays containing the start values
    for box in [lolat, hilat, deep]:
        for k in model_vars:
            box[k] = np.full(time.shape, box[k])
    print(lolat)
    fluxes = {}  # Create a dictionary to keep track of the fluxes calculated at each step

    ### LOOP STARTS HERE ###
    for i in range(1, time.size):
        last = i - 1  # get the index of the previous time step

        # 1. calculate the thermohaline circulation flux, Q_T

        Q_T = Q_k*(Q_alpha*(lolat['T'][last] - hilat['T'][last]) - Q_beta*(lolat['S'][last] - hilat['S'][last]))
        fluxes['Q_T'] = Q_T
        print(Q_T)

        # 2. calculate all the fluxes
        # Note: be careful with the timestep parameter, dt - how do you include this in the fluxes?

        # 2.i calculate the mixing fluxes for each model variable
        # for var in model_vars:
        #   calculate each of the components of the mixing flux for each variable in turn, storing the output in the fluxes dictionary
        #   for example, for the thermohaline circulation flux, you could write:
        #   fluxes[f'Q_{var}_deep'] = Q_T * (hilat[var][last] - deep[var][last]) * dt  # mol dt-1

        #   [[[do it here...]]]

        # 2.ii calculate temperature exchange with each surface box
        # for box in [lolat, hilat]:
        #   calculate the temperature exchange flux for each surface box, storing the output in the fluxes dictionary

        #   [[[do it here...]]]

        # 3. use the calculated fluxes to update the state of the model variables in each box

        # 3.i apply fluxes to calculate new values in deep box
        # for var in model_vars:
        #   apply the calculated fluxes to update the state of each model variable in the deep box

        #   [[[do it here...]]]

        # 3.ii apply fluxes to calculate new values in surface boxes
        # for box in [lolat, hilat]:
        #   apply the calculated fluxes to update the state of each model variable in each surface box
        #   (N.B. you can't loop through the model variables here, as the fluxes are different for each box)

        #   [[[do it here...]]]

    return time, lolat, hilat, deep

ocean_model(init_lolat, init_hilat, init_deep, 1000, 0.5)