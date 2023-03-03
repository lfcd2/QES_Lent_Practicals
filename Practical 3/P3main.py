import numpy as np
import matplotlib.pyplot as plt
from OceanTools.tools import plot

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
    'E': E  # salt added due to evaporation - precipitation, kg m-3 yr-1
}
init_lolat['V'] = init_lolat['SA'] *  init_lolat['depth']  # box volume, m3

init_deep = {
    'name': 'deep',
    'SA': SA_ocean,  # box surface area, m2
    'T': 5.,  # initial water temperature, Celcius
    'S': 34.5,  # initial salinity
    'V': V_ocean - init_lolat['V'] - init_hilat['V']
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


        # 2. calculate all the fluxes
        # Note: be careful with the timestep parameter, dt - how do you include this in the fluxes?
        for var in model_vars:
            fluxes[f'Q_{var}_deep'] = Q_T * (hilat[var][last] - deep[var][last]) * dt  # mol dt-1
            fluxes[f'Q_{var}_hilat'] = Q_T * (lolat[var][last] - hilat[var][last]) * dt  # mol dt-1
            fluxes[f'Q_{var}_lolat'] = Q_T * (deep[var][last] - lolat[var][last]) * dt  # mol dt-1

            fluxes[f'vmix_{var}_hilat'] = hilat['V'] / hilat['tau_M'] * (hilat[var][last] - deep[var][last]) * dt  # mol dt-1
            fluxes[f'vmix_{var}_lolat'] = lolat['V'] / lolat['tau_M'] * (lolat[var][last] - deep[var][last]) * dt  # mol dt-1

        # 2.i calculate the mixing fluxes for each model variable
        for box in [hilat, lolat]:
            boxname = box['name']
            fluxes[f'dT_{boxname}'] = box['V'] / box['tau_T'] * (box['T_atmos'] - box['T'][last]) * dt  # mol dt-1

        #   [[[do it here...]]]

        # 2.ii calculate temperature exchange with each surface box
        for var in model_vars:
            deep[var][i] = deep[var][last] + (
                    fluxes[f'Q_{var}_deep'] + fluxes[f'vmix_{var}_hilat'] + fluxes[f'vmix_{var}_lolat']) / deep['V']

        # 3. use the calculated fluxes to update the state of the model variables in each box
        for box in [hilat, lolat]:
            boxname = box['name']
            box['S'][i] = box['S'][last] + (fluxes[f'Q_S_{boxname}'] - fluxes[f'vmix_S_{boxname}'] + box['E'] * dt) / box['V']
            box['T'][i] = box['T'][last] + (fluxes[f'Q_T_{boxname}'] - fluxes[f'vmix_T_{boxname}'] + fluxes[f'dT_{boxname}']) / box['V']

    return time, lolat, hilat, deep

time, lolat, hilat, deep = ocean_model(init_lolat, init_hilat, init_deep, 1000, 0.5)

fig, axs = plot.boxes(time, ['T', 'S'], lolat, hilat, deep)
plt.show()