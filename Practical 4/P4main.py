import cbsyst
import numpy as np
import matplotlib.pyplot as plt
from cbsyst import Csys

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

# set up boxes
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
    'DIC': (38700e15/12)/3, # TODO is this correct units?
    'TA': 3.1e18,
    'tau_CO2': 2.
}
init_hilat['V'] = init_hilat['SA'] *  init_hilat['depth']  # box volume, m3

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
    # Add new variables here
    'DIC': (38700e15/12)/3, # TODO is this correct units?
    'TA': 1.,
    'tau_CO2': 2.
}
init_lolat['V'] = init_lolat['SA'] *  init_lolat['depth']  # box volume, m3

init_deep = {
    'name': 'deep',
    'V': V_ocean - init_lolat['V'] - init_hilat['V'],  # box volume, m3
    'T': 5.483637,  # initial water temperature, Celcius
    'S': 34.47283,  # initial salinity
    # Add new variables here
    'DIC': (38700e15/12)/3, # TODO is this correct units?
    'TA': 1.,
}

# create a new dictionary for the atmosphere here, then calculate the total moles of CO2 in the atmosphere
init_atmos = {
    'name': 'atmos',
    'mass': 5e21, # in grams
    'moles_air': 1.736e20,
    'moles_CO2': 850e15/12,
    'GtC_emissions': 0
}
init_atmos['pCO2'] = 1e6*init_atmos['moles_CO2']/init_atmos['moles_air']


def ocean_model(lolat, hilat, deep, atmos, tmax, dt):
    """Run the ocean model for a given time period and return the results for each box.

    Parameters
    ----------
    lolat, hilat, deep, atmos : dict
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
    model_vars = ['T', 'S', 'DIC', 'TA']
    atmos_model_vars = ['moles_CO2', 'pCO2']

    # create copies of the input dictionaries so we don't modify the originals
    lolat = lolat.copy()
    hilat = hilat.copy()
    deep = deep.copy()
    atmos = atmos.copy()

    # turn all time-evolving variables into arrays containing the start values
    for box in [lolat, hilat, deep]:
        for k in model_vars:
            box[k] = np.full(time.shape, box[k])
    for k in atmos_model_vars:
        atmos[k] = np.full(time.shape, atmos[k])
    if isinstance(atmos['GtC_emissions'], (int, float)):
        atmos['GtC_emissions'] = np.full(time.shape, atmos['GtC_emissions'])

    # calculate initial surface carbon chemistry in the surface boxes using Csys, and store a few key variables - CO2, pH, pCO2 and K0
    # NOTE: you'll want to re-use this code lower down to re-calculated carbon speciation after you've updated DIC and TA in the surface.
    for box in [lolat, hilat]:
        csys = Csys(
            TA=1e3 * box['TA'],  # 1e3 is necessary because cbsyst expects units of umol L-1
            DIC=1e3 * box['DIC'],  # 1e3 is necessary because cbsyst expects units of umol L-1
            T_in=box['T'], S_in=box['S'],
        )
        box['CO2'] = csys.CO2 * 1e-3  # 1e3 converts back to mol m-3
        box['pH'] = csys.pHtot  # this is -log10([H+])
        box['pCO2'] = csys.pCO2  # this is already in the right units (ppm)
        box['K0'] = csys.Ks.K0  # Henry's law constant - [CO2] = K0 * pCO2



    # Create a dictionary to keep track of the fluxes calculated at each step
    fluxes = {}

    ### LOOP STARTS HERE ###
    for i in range(1, time.size):
        last = i - 1  # index of last model step

        # calculate circulation flux, Q
        dT = lolat['T'][last] - hilat['T'][last]
        dS = lolat['S'][last] - hilat['S'][last]
        Q = Q_k * (Q_alpha * dT - Q_beta * dS)

        # calculate mixing fluxes for model variables
        for var in model_vars:
            # Nothing to do here! If you've added DIC to the model_vars list above, the mixing and circulation of DIC are calculated automatically here
            fluxes[f'Q_{var}_deep'] = Q * (hilat[var][last] - deep[var][last]) * dt  # amount dt-1
            fluxes[f'Q_{var}_hilat'] = Q * (lolat[var][last] - hilat[var][last]) * dt  # amount dt-1
            fluxes[f'Q_{var}_lolat'] = Q * (deep[var][last] - lolat[var][last]) * dt  # amount dt-1

            fluxes[f'vmix_{var}_hilat'] = hilat['V'] / hilat['tau_M'] * (
                        hilat[var][last] - deep[var][last]) * dt  # amount dt-1
            fluxes[f'vmix_{var}_lolat'] = lolat['V'] / lolat['tau_M'] * (
                        lolat[var][last] - deep[var][last]) * dt  # amount dt-1

        # calculate surface-specific fluxes
        for box in [hilat, lolat]:
            boxname = box['name']
            # temperature exchange with atmosphere
            fluxes[f'dT_{boxname}'] = box['V'] / box['tau_T'] * (box['T_atmos'] - box['T'][last]) * dt  # mol dt-1

            # TODO: Calculate the fluxes for CO2 between the atmosphere and each surface box here
            # NOTE: be careful with units! Your CO2 change should be in units of mol m-3.
            fluxes[f'dDIC_{boxname}'] = (box['V'] / box['tau_CO2']) * (box['DIC'][last] - box['K0'][last]*box['pCO2'][last]) * dt

        # TODO: calculate the flux of CO2 generated by emissions (this will be zero to start with, because GtC_emissions is zero above, but it will become useful later on when you simulate adding CO2 into the atmosphere)
        fluxes['emissions'] = atmos['GtC_emissions']*1e15/12

        # update deep box
        for var in model_vars:
            deep[var][i] = (deep[var][last] + (fluxes[f'Q_{var}_deep'] + fluxes[f'vmix_{var}_hilat'] + fluxes[f'vmix_{var}_lolat']) / deep['V'])
            # Nothing to do here! You've added 'DIC' to the 'model_vars' list above, so it is already updated here.

        # update surface boxes
        for box in [hilat, lolat]:
            boxname = box['name']
            box['S'][i] = box['S'][last] + (fluxes[f'Q_S_{boxname}'] - fluxes[f'vmix_S_{boxname}'] + box['E'] * dt) / box['V']  # salinity dt-1
            box['T'][i] = box['T'][last] + (fluxes[f'Q_T_{boxname}'] - fluxes[f'vmix_T_{boxname}'] + fluxes[f'dT_{boxname}']) / box['V']  # degrees dt-1
            # --> TODO: Add DIC and TA here

            box['DIC'][i] = box['DIC'][last] + (fluxes[f'Q_DIC_{boxname}'] - fluxes[f'vmix_DIC_{boxname}'] + fluxes[f'dDIC_{boxname}']) / box['V']

            print(box['DIC'][last], fluxes[f'Q_DIC_{boxname}'], fluxes[f'vmix_DIC_{boxname}'], fluxes[f'dDIC_{boxname}'])
            box['TA'][i] = box['TA'][last] + (fluxes[f'Q_TA_{boxname}'] - fluxes[f'vmix_TA_{boxname}']) / box['V']
            # TODO: calculate carbon speciation in the surface boxes (use the example code above)

        # TODO: update CO2 in the atmosphere
        # NOTE: be careful with units, and calculate both moles of CO2 in the atmosphere and the pCO2.
        atmos['moles_CO2'][i] = atmos['moles_CO2'][last] + fluxes['emissions'][last]
        atmos['pCO2'][i] = 1e6 * atmos['moles_CO2'][i] / init_atmos['moles_air']

    return time, lolat, hilat, deep, atmos

time, lolat, hilat, deep, atmos = ocean_model(init_lolat, init_hilat, init_deep, init_atmos, 3000, 0.5)
from OceanTools.tools import plot

fig, axs = plot.boxes(time, ['DIC', 'TA', 'pCO2'], lolat, hilat, deep, atmos)
print(lolat['DIC'])
plt.show()