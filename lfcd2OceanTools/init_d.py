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
