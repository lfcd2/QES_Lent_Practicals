import numpy as np
import matplotlib.pyplot as plt
from cbsyst import Csys
from tqdm import tqdm  # This isn't strictly speaking necessary, but it makes it easier to track runtime progress
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
E = Fw * SA_ocean * (1 - fSA_hilat) * Sref  # amount of salt removed from the low lat box, g kg-1 yr-1, ~ kg m-3 yr-1


def initialise_dicts_15():
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
        'tau_CO2': 2.,  # timescale of CO2 exchange, yr
        'DIC': 2.32226,  # Dissolved Inorganic Carbon concentration, mol m-3
        'TA': 3.1e18/V_ocean,  # Total Alkalinity, mol m-3
        'tau_PO4': 3.,
        'f_CaCO3': 0.2,
        'PO4': 3e15/V_ocean
    }
    init_hilat['V'] = init_hilat['SA'] * init_hilat['depth']  # box volume, m3

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
        'tau_CO2': 2.,  # timescale of CO2 exchange, yr
        'DIC': 2.26201,  # Dissolved Inorganic Carbon concentration, mol m-3
        'TA': 3.1e18/V_ocean,  # Total Alkalinity, mol m-3
        'tau_PO4': 2.,
        'f_CaCO3': 0.3,
        'PO4': 3e15 / V_ocean
    }
    init_lolat['V'] = init_lolat['SA'] * init_lolat['depth']  # box volume, m3

    init_deep = {
        'name': 'deep',
        'V': V_ocean - init_lolat['V'] - init_hilat['V'],  # box volume, m3
        'T': 5.483637,  # initial water temperature, Celcius
        'S': 34.47283,  # initial salinity
        'DIC': 2.32207,  # Dissolved Inorganic Carbon concentration, mol m-3
        'TA': 3.1e18/V_ocean,  # Total Alkalinity, mol m-3
        'PO4': 3e15/V_ocean
    }

    init_atmos = {
        'name': 'atmos',
        'mass': 5.132e18,  # kg
        'moles_air': 1.736e20,  # moles
        'moles_CO2': 850e15 / 12,  # moles
        'GtC_emissions': 0.0  # annual emissions of CO2 into the atmosphere, GtC
    }
    init_atmos['pCO2'] = init_atmos['moles_CO2'] / init_atmos['moles_air'] * 1e6

    # returns the dictionaries
    return init_lolat, init_hilat, init_deep, init_atmos


def ocean_model_p15(dicts, tmax, dt):
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

            # Productivity
            fluxes[f'prod_PO4_{box_name}'] = box['V'] / box['tau_PO4'] * box['PO4'][last] * dt

            # productivity on DIC (we add 106 here because the bio pump is favourable for DIC, but unfavourable for TA)
            fluxes[f'prod_DIC_{box_name}'] = fluxes[f'prod_PO4_{box_name}'] * (106 * box['f_CaCO3'] + 106)

            # productivity on TA
            fluxes[f'prod_TA_{box_name}'] = fluxes[f'prod_PO4_{box_name}'] * (106 * box['f_CaCO3'] * 2 - 18)

            # calculate emissions flux
        fluxes['emissions'] = atmos['GtC_emissions'] * 1e15 / 12 * dt

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
                TA=1e3 * box['TA'][i],  # 1e3 is necessary because cbsyst expects units of umol L-1
                DIC=1e3 * box['DIC'][i],  # 1e3 is necessary because cbsyst expects units of umol L-1
                T_in=box['T'][i], S_in=box['S'][i],
            )
            box['CO2'][i] = csys.CO2 * 1e-3  # 1e3 converts back to mol m-3
            box['pH'][i] = csys.pHtot  # this is -log10([H+])
            box['pCO2'][i] = csys.pCO2  # this is already in the right units (ppm)
            box['K0'][i] = csys.Ks.K0  # Henry's law constant - [CO2] = K0 * pCO2

        # ===================== UPDATE ATMOSPHERE BOX ===================== #

        # NOTE: be careful with units, and calculate both moles of CO2 in the atmosphere and the pCO2.
        atmos['moles_CO2'][i] = (atmos['moles_CO2'][last]
                                 + fluxes[f'dDIC_hilat']
                                 + fluxes[f'dDIC_lolat']
                                 + fluxes['emissions'][last])

        atmos['pCO2'][i] = 1e6 * atmos['moles_CO2'][i] / atmos['moles_air']

    '''
    This is the end of the loop
    ###################################################################################################################
    '''

    return time, (lolat, hilat, deep, atmos)


def run():
    """
    this runs the code
    """

    # this line of code runs initialise_dicts which sets up the initial dictionaries
    dicts = initialise_dicts_15()

    # this line of code runs the model
    time_array, finished_dicts = ocean_model_p15(dicts, 1000, 0.5)

    # this unpacks the result that is output from the model
    final_lolat, final_hilat, final_deep, final_atmos = finished_dicts

    # this uses oscars plot function to plot DIC, TA and pCO2 (you can experiment by adding variables into this list)
    fig, axs = plot.boxes(time_array, ['DIC', 'TA', 'pCO2', 'PO4'],
                          final_lolat, final_hilat, final_deep, final_atmos)

    # adjust axes of the pCO2 plot
    axs[-2].set_ylim(0, 1600)

    # plot the graph
    plt.show()


if __name__ == '__main__':
    run()
