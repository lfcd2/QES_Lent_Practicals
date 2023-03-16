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
rho_org = 1100
rho_CaCO3 = 2700


def ocean_model_q3(dicts, tmax, dt):
    """Run the ocean model for a given time period and return the results for each box.

    Parameters
    ----------
    dicts : list
        list of dictionaries containing the box properties
    tmax : int or float
        The maximum time to run the model for (yr)
    dt : float
        The time step of the model (yr)

    Returns
    -------
    time : list of time steps
    tuple of (lolat, hilat, deep, atmos)
    """

    # unpack the dicts from the list to reduce the number of variables we input into the function
    lolat, hilat, deep, atmos = dicts

    # create the timescale for the model - we will have to multiply all the fluxes by this value
    time = np.arange(0, tmax + dt, dt)

    # identify which variables will change with time - we will iterate over these lists
    model_vars = ['T', 'S', 'DIC', 'TA', 'PO4']
    atmos_model_vars = ['moles_CO2', 'pCO2']
    track_vars = ['f_CaCO3', 'particle_sinking_time', 'v']

    '''
    create copies of the input dictionaries, so we don't modify the originals
    this is unnecessary if you are using a decent IDE and running the code locally, so I've commented them out
    '''
    # lolat = lolat.copy()
    # hilat = hilat.copy()
    # deep = deep.copy()
    # atmos = atmos.copy()

    '''
    turn all time-evolving variables into arrays containing the start values
    this iterates over each box and variable, converting it from a single number into a list of that number with same
    lengths as the number of time steps you are doing your model for - it does it first for the model_vars and then
    for the atmos_model_vars
    '''
    for box in [lolat, hilat, deep]:
        for k in model_vars:
            box[k] = np.full(time.shape, box[k])
        for k in track_vars:
            if k in box:
                box[k] = np.full(time.shape, box[k])

    for k in atmos_model_vars:
        atmos[k] = np.full(time.shape, atmos[k])
    if isinstance(atmos['GtC_emissions'], (int, float)):
        atmos['GtC_emissions'] = np.full(time.shape, atmos['GtC_emissions'])

    '''
    Uses Oscars cbsyst package to calculate initial surface carbon chemistry in the surface boxes using Csys, and store
    a few key variables: CO2, pH, pCO2 and K0
    '''
    for box in [lolat, hilat]:
        csys = Csys(
            TA=box['TA'],
            DIC=box['DIC'],
            T_in=box['T'], S_in=box['S'],
            unit='mmol'
        )
        box['CO2'] = csys.CO2
        box['pH'] = csys.pHtot  # this is -log10([H+])
        box['pCO2'] = csys.pCO2  # this is already in the right units (ppm)
        box['K0'] = csys.Ks.K0  # Henry's law constant - [CO2] = K0 * pCO2
        box['CO3'] = csys.CO3
        box['Omega'] = csys.OmegaA

        f_remaining = np.exp(k_diss * box['particle_sinking_time'][0] * (Omega_crit - box['Omega']) ** n_diss)
        f_remaining[box['Omega'] > Omega_crit] = 1
        box['f_CaCO3'] = calc_slope * box['Omega'] * f_remaining

    # Create a dictionary to keep track of the fluxes calculated at each step
    fluxes = {}

    '''
    ###################################################################################################################
    This is where the loop starts - this repeats for each timestep. It calculates the fluxes, adds them to the previous
    values, then uses oscars package to recalculate DIC, pCO2 and the variables you need to calculate TA
    '''
    for i in tqdm(range(1, time.size)):
        last = i - 1  # index of last model step

        # ===================== CALCULATE FLUXES ===================== #

        # calculate circulation flux, Q
        dT = lolat['T'][last] - hilat['T'][last]
        dS = lolat['S'][last] - hilat['S'][last]
        Q = Q_k * (Q_alpha * dT - Q_beta * dS)

        # calculate mixing fluxes for model variables, all in units of amount dt-1
        for var in model_vars:
            fluxes[f'Q_{var}_deep'] = Q * (hilat[var][last] - deep[var][last]) * dt
            fluxes[f'Q_{var}_hilat'] = Q * (lolat[var][last] - hilat[var][last]) * dt
            fluxes[f'Q_{var}_lolat'] = Q * (deep[var][last] - lolat[var][last]) * dt

            fluxes[f'vmix_{var}_hilat'] = hilat['V'] / hilat['tau_M'] * (hilat[var][last] - deep[var][last]) * dt
            fluxes[f'vmix_{var}_lolat'] = lolat['V'] / lolat['tau_M'] * (lolat[var][last] - deep[var][last]) * dt

        # calculate surface-specific fluxes
        for box in [hilat, lolat]:
            box_name = box['name']

            # temperature exchange with atmosphere
            fluxes[f'dT_{box_name}'] = \
                box['V'] / box['tau_T'] * (box['T_atmos'] - box['T'][last]) * dt  # mol dt-1

            # DIC exchange with atmosphere (careful with units! Your CO2 change should be in units of mol m-3.)
            fluxes[f'dDIC_{box_name}'] = \
                (box['V'] / box['tau_CO2']) * (box['CO2'][last] - box['K0'][last] * atmos['pCO2'][last] * 1e-3) * dt

            rho_p = (rho_org + box['f_CaCO3'][last] * 10 / 3 * rho_org) / (1 + box['f_CaCO3'][last] * 10 / 3 * rho_org / rho_CaCO3)
            v = particle_velocity * (rho_p - 1000) / (box['rho_particle'] - 1000)
            box['particle_sinking_time'][i] = box['depth'] / v
            exponential_factor = np.exp(- box['k_ballast'] * box['particle_sinking_time'][i] / (24*60*60))

            # Productivity
            fluxes[f'prod_PO4_{box_name}'] = box['V'] / box['tau_PO4'] * box['PO4'][last] * dt * 5 * exponential_factor

            # productivity on DIC (we add 106 here because the bio pump is favourable for DIC, but unfavourable for TA)
            fluxes[f'prod_DIC_{box_name}'] = fluxes[f'prod_PO4_{box_name}'] * (106 * box['f_CaCO3'][last] + 106)

            # productivity on TA
            fluxes[f'prod_TA_{box_name}'] = fluxes[f'prod_PO4_{box_name}'] * (106 * box['f_CaCO3'][last] * 2 - 18)

            # calculate emissions flux
        fluxes['emissions'] = atmos['GtC_emissions'][last] * 1e15 / 12 * dt

        # ===================== UPDATE BOXES ===================== #

        # update deep box
        for var in model_vars:
            deep[var][i] = deep[var][last] + (fluxes[f'Q_{var}_deep']
                                              + fluxes[f'vmix_{var}_hilat']
                                              + fluxes[f'vmix_{var}_lolat']) / deep['V']
            if var in ['DIC', 'TA', 'PO4']:  # this takes into account  the bio pump for the relevant vars
                for box_loc in ['hilat', 'lolat']:
                    deep[var][i] += fluxes[f'prod_{var}_{box_loc}'] / deep['V']

        # update surface boxes for each variable
        for box in [hilat, lolat]:
            box_name = box['name']
            box['S'][i] = box['S'][last] + (fluxes[f'Q_S_{box_name}']
                                            - fluxes[f'vmix_S_{box_name}']
                                            + box['E'] * dt) / box['V']  # salinity dt-1

            box['T'][i] = box['T'][last] + (fluxes[f'Q_T_{box_name}']
                                            - fluxes[f'vmix_T_{box_name}']
                                            + fluxes[f'dT_{box_name}']) / box['V']  # degrees dt-1

            box['DIC'][i] = box['DIC'][last] + (fluxes[f'Q_DIC_{box_name}']
                                                - fluxes[f'vmix_DIC_{box_name}']
                                                - fluxes[f'dDIC_{box_name}']
                                                - fluxes[f'prod_DIC_{box_name}']) / box['V']  # mol m-3 dt-1

            box['TA'][i] = box['TA'][last] + (fluxes[f'Q_TA_{box_name}']
                                              - fluxes[f'vmix_TA_{box_name}']
                                              - fluxes[f'prod_TA_{box_name}']) / box['V']  # mol m-3 dt-1

            box['PO4'][i] = box['PO4'][last] + (fluxes[f'Q_PO4_{box_name}']
                                                - fluxes[f'vmix_PO4_{box_name}']
                                                - fluxes[f'prod_PO4_{box_name}']) / box['V']  # mol m-3 dt-1

        # ===================== RECALCULATE THE DIC ETC FOR EACH BOX ===================== #

        for box in [lolat, hilat]:
            csys = Csys(
                TA=box['TA'][i],
                DIC=box['DIC'][i],
                T_in=box['T'][i], S_in=box['S'][i],
                unit='mmol'
            )
            box['CO2'][i] = csys.CO2
            box['pH'][i] = csys.pHtot  # this is -log10([H+])
            box['pCO2'][i] = csys.pCO2  # this is already in the right units (ppm)
            box['K0'][i] = csys.Ks.K0  # Henry's law constant - [CO2] = K0 * pCO2
            box['CO3'][i] = csys.CO3
            box['Omega'][i] = csys.OmegaA

            # update f_CaCO3
            if box['Omega'][i] > Omega_crit:
                f_remaining = 1
            else:
                f_remaining = np.exp(k_diss * box['particle_sinking_time'][last] * (Omega_crit - box['Omega'][i]) ** n_diss)
            box['f_CaCO3'][i] = calc_slope * box['Omega'][i] * f_remaining

        # ===================== UPDATE ATMOSPHERE BOX ===================== #

        # NOTE: be careful with units, and calculate both moles of CO2 in the atmosphere and the pCO2.
        atmos['moles_CO2'][i] = (atmos['moles_CO2'][last]
                                 + fluxes[f'dDIC_hilat']
                                 + fluxes[f'dDIC_lolat']
                                 + fluxes['emissions'])

        atmos['pCO2'][i] = 1e6 * atmos['moles_CO2'][i] / atmos['moles_air']

    '''
    This is the end of the loop
    ###################################################################################################################
    '''

    return time, (lolat, hilat, deep, atmos)
