from cbsyst import Csys
import numpy as np
import matplotlib.pyplot as plt
from OceanTools.tools import plot
from OceanTools.tools.helpers import get_last_values
from tqdm import tqdm

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
E = Fw * SA_ocean * (1 - fSA_hilat) * Sref
# amount of salt removed from the low latitude box,  g kg-1 yr-1, ~ kg m-3 yr-1

particle_velocity = 10  # m d-1
k_diss = -0.07  # d-1
n_diss = 2.0  # unitless
Omega_crit = 2.5  # unitless
calc_slope = 0.12  # f_CaCO3 / Omega


def acidification_model(dicts, tmax, dt):
    lolat, hilat, deep, atmos = dicts
    # create the time scale for the model
    time = np.arange(0, tmax + dt, dt)

    # identify which variables will change with time
    model_vars = ['T', 'S', 'DIC', 'TA', 'PO4']
    atmos_model_vars = ['moles_CO2', 'pCO2']
    ### NEW CODE
    track_vars = ['f_CaCO3', 'particle_sinking_time', 'v', 'DIC_export', 'TA_export']
    ###

    # create copies of the input dictionaries so we don't modify the originals
    lolat = lolat.copy()
    hilat = hilat.copy()
    deep = deep.copy()
    atmos = atmos.copy()

    # turn all time-evolving variables into arrays containing the start values
    for box in [lolat, hilat, deep]:
        if box != deep:
            box['DIC_export'] = 0.
            box['TA_export'] = 0.
        for k in model_vars:
            box[k] = np.full(time.shape, box[k])
        ### NEW CODE
        for k in track_vars:
            if k in box:
                box[k] = np.full(time.shape, box[k])
        ###
    for k in atmos_model_vars:
        atmos[k] = np.full(time.shape, atmos[k])
    if isinstance(atmos['GtC_emissions'], (int, float)):
        atmos['GtC_emissions'] = np.full(time.shape, atmos['GtC_emissions'])

    # calculate initial surface carbon chemistry in the surface boxes using Csys, and store a few key variables - CO2, pH, pCO2 and K0
    for box in [lolat, hilat]:
        csys = Csys(
            TA=box['TA'],
            DIC=box['DIC'],
            T_in=box['T'], S_in=box['S'],
            unit='mmol'
        )
        box['CO2'] = csys.CO2
        box['pH'] = csys.pHtot
        box['pCO2'] = csys.pCO2
        box['K0'] = csys.Ks.K0

        ### NEW CODE
        box['CO3'] = csys.CO3
        box['Omega'] = csys.OmegaA

        # calculate initial f_CaCO3
        f_remaining = np.exp(k_diss * box['particle_sinking_time'][0] * (Omega_crit - box['Omega']) ** n_diss)
        f_remaining[box['Omega'] > Omega_crit] = 1
        box['f_CaCO3'] = calc_slope * box['Omega'] * f_remaining
        ###

    # Create a dictionary to keep track of the fluxes calculated at each step
    fluxes = {}

    for i in range(1, time.size):
        last = i - 1  # index of last model step

        # calculate circulation flux, Q
        dT = lolat['T'][last] - hilat['T'][last]
        dS = lolat['S'][last] - hilat['S'][last]
        Q = Q_k * (Q_alpha * dT - Q_beta * dS)

        # calculate mixing fluxes for model variables
        for var in model_vars:
            fluxes[f'Q_{var}_deep'] = Q * (hilat[var][last] - deep[var][last]) * dt  # mol dt-1
            fluxes[f'Q_{var}_hilat'] = Q * (lolat[var][last] - hilat[var][last]) * dt  # mol dt-1
            fluxes[f'Q_{var}_lolat'] = Q * (deep[var][last] - lolat[var][last]) * dt  # mol dt-1

            fluxes[f'vmix_{var}_hilat'] = hilat['V'] / hilat['tau_M'] * (
                        hilat[var][last] - deep[var][last]) * dt  # mol dt-1
            fluxes[f'vmix_{var}_lolat'] = lolat['V'] / lolat['tau_M'] * (
                        lolat[var][last] - deep[var][last]) * dt  # mol dt-1

        # calculate surface-specific fluxes
        for box in [hilat, lolat]:
            boxname = box['name']

            # temperature exchange with atmosphere
            fluxes[f'dT_{boxname}'] = box['V'] / box['tau_T'] * (box['T_atmos'] - box['T'][last]) * dt  # mol dt-1
            # CO2 exchange with atmosphere
            fluxes[f'dCO2_{boxname}'] = box['V'] / box['tau_CO2'] * (
                        box['CO2'][last] - 1e-3 * atmos['pCO2'][last] * box['K0'][last]) * dt  # mol dt-1

            # organic matter production
            fluxes[f'export_PO4_{boxname}'] = box['PO4'][last] * box['V'] / box['tau_PO4'] * dt # mol PO4 dt-1

            ### MODIFIED CODE
            # DIC export by productivity :                                  redfield + calcification
            fluxes[f'export_DIC_{boxname}'] = fluxes[f'export_PO4_{boxname}'] * (
                        106 + 106 * box['f_CaCO3'][last])  # mol DIC dt-1
            box['DIC_export'][i] = fluxes[f'export_DIC_{boxname}'] / (dt*box['V'])

            # TA export by productivity :                                  redfield + calcification
            fluxes[f'export_TA_{boxname}'] = fluxes[f'export_PO4_{boxname}'] * (
                        -18 + 2 * 106 * box['f_CaCO3'][last])  # mol TA dt-1
            box['TA_export'][i] = fluxes[f'export_TA_{boxname}'] / (dt*box['V'])

        fluxes['dCO2_emissions'] = atmos['GtC_emissions'][last] * 1e15 / 12 * dt  # mol dt-1

        # update deep box
        for var in model_vars:
            if var in ['T', 'S']:
                deep[var][i] = deep[var][last] + (
                        fluxes[f'Q_{var}_deep'] + fluxes[f'vmix_{var}_hilat'] + fluxes[f'vmix_{var}_lolat']
                ) / deep['V']
            else:
                deep[var][i] = deep[var][last] + (
                        fluxes[f'Q_{var}_deep'] + fluxes[f'vmix_{var}_hilat'] + fluxes[f'vmix_{var}_lolat'] + fluxes[
                    f'export_{var}_hilat'] + fluxes[f'export_{var}_lolat']
                ) / deep['V']

        # update surface boxes
        for box in [hilat, lolat]:
            boxname = box['name']
            box['S'][i] = box['S'][last] + (fluxes[f'Q_S_{boxname}'] - fluxes[f'vmix_S_{boxname}'] + box['E'] * dt) / \
                          box['V']
            box['T'][i] = box['T'][last] + (
                        fluxes[f'Q_T_{boxname}'] - fluxes[f'vmix_T_{boxname}'] + fluxes[f'dT_{boxname}']) / box['V']
            box['DIC'][i] = box['DIC'][last] + (
                        fluxes[f'Q_DIC_{boxname}'] - fluxes[f'vmix_DIC_{boxname}'] - fluxes[f'dCO2_{boxname}'] - fluxes[
                    f'export_DIC_{boxname}']) / box['V']
            box['TA'][i] = box['TA'][last] + (
                        fluxes[f'Q_TA_{boxname}'] - fluxes[f'vmix_TA_{boxname}'] - fluxes[f'export_TA_{boxname}']) / \
                           box['V']
            box['PO4'][i] = box['PO4'][last] + (
                        fluxes[f'Q_PO4_{boxname}'] - fluxes[f'vmix_PO4_{boxname}'] - fluxes[f'export_PO4_{boxname}']) / \
                            box['V']

            # update carbon speciation
            csys = Csys(
                TA=box['TA'][i],
                DIC=box['DIC'][i],
                T_in=box['T'][i], S_in=box['S'][i],
                unit='mmol'
            )
            box['CO2'][i] = csys.CO2
            box['pCO2'][i] = csys.pCO2
            box['pH'][i] = csys.pHtot
            box['K0'][i] = csys.Ks.K0

            ### NEW CODE
            box['CO3'][i] = csys.CO3
            box['Omega'][i] = csys.OmegaA

            # update f_CaCO3
            if box['Omega'][i] > Omega_crit:
                f_remaining = 1
            else:
                f_remaining = np.exp(k_diss * box['particle_sinking_time'][0] * (Omega_crit - box['Omega'][i]) ** n_diss)
            box['f_CaCO3'][i] = calc_slope * box['Omega'][i] * f_remaining
            ###

        # update atmosphere
        atmos['moles_CO2'][i] = atmos['moles_CO2'][last] + fluxes['dCO2_hilat'] + fluxes['dCO2_lolat'] + fluxes[
            'dCO2_emissions']
        atmos['pCO2'][i] = 1e6 * atmos['moles_CO2'][i] / atmos['moles_air']

    return time, (lolat, hilat, deep, atmos)