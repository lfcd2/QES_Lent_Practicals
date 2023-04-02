import numpy as np
import matplotlib.pyplot as plt
from OceanTools.tools import plot
from newOceanTools.newTools import copy_dicts, add_emissions, add_fancy_labels

from Model1Original import original_model
from Model2Acidification import acidification_model
from Model3Ballasting import ballasting_model


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
E = Fw * SA_ocean * (1 - fSA_hilat) * Sref  # amount of salt removed from the low lat box,  g kg-1 yr-1, ~ kg m-3 yr-1

particle_velocity = 10  # m d-1
k_diss = -0.07  # d-1
n_diss = 2.0  # unitless
Omega_crit = 2.5  # unitless
calc_slope = 0.12  # f_CaCO3 / Omega
rho_org = 1200
rho_CaCO3 = 2700


def init_dicts_16():
    # NOTE: Initial DIC, TA, PO4 and pCO2 values are set to steady state values from the Ocean Acidification model.

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
        'DIC': 2.02823,  # Dissolved Inorganic Carbon concentration, mol m-3
        'TA': 2.22043,  # Total Alkalinity, mol m-3
        'tau_PO4': 3.,  # phosphate half life, yr at initial f_CaCO3
        'PO4': 8.90099e-05,  # Phosphate conc, mol m-3
        'f_CaCO3': 0.18134,  # fraction of organic matter export that produces CaCO3 at starting [CO3]
        'k_ballast': 0.0805,
        'rho_particle': 1416.27,
        'exp': 0.2
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
        'DIC': 1.99301,  # Dissolved Inorganic Carbon concentration, mol m-3
        'TA': 2.21683,  # Total Alkalinity, mol m-3
        'tau_PO4': 2.,  # phosphate half life, yr at initial f_CaCO3
        'PO4': 1.65460e-04,  # Phosphate conc, mol m-3
        'f_CaCO3': 0.30453,  # fraction of organic matter export that produces CaCO3 at starting [CO3]
        'k_ballast': 0.1609,
        'rho_particle': 1568.14,
        'exp': 0.2
    }
    init_lolat['V'] = init_lolat['SA'] * init_lolat['depth']  # box volume, m3

    init_deep = {
        'name': 'deep',
        'V': V_ocean - init_lolat['V'] - init_hilat['V'],  # box volume, m3
        'T': 5.483637,  # initial water temperature, Celcius
        'S': 34.47283,  # initial salinity
        'DIC': 2.32710,  # Dissolved Inorganic Carbon concentration, mol m-3
        'TA': 2.31645,  # Total Alkalinity, mol m-3
        'PO4': 2.30515e-03,  # Phosphate conc, mol m-3
    }

    init_atmos = {
        'name': 'atmos',
        'mass': 5.132e18,  # kg
        'moles_air': 1.736e20,  # moles
        'moles_CO2': 872e15 / 12,  # moles
        'GtC_emissions': 0.0  # annual emissions of CO2 into the atmosphere, GtC
    }
    init_atmos['pCO2'] = init_atmos['moles_CO2'] / init_atmos['moles_air'] * 1e6

    init_hilat['particle_sinking_time'] = init_hilat['depth'] / particle_velocity
    init_lolat['particle_sinking_time'] = init_lolat['depth'] / particle_velocity

    return [init_lolat, init_hilat, init_deep, init_atmos]


def plot_all_three():
    """
    this runs the code
    """

    # this line of code runs initialise_dicts which sets up the initial dictionaries
    init_dicts = init_dicts_16()
    tmax = 2000  # how many years to simulate (yr)
    dt = 0.5  # the time step of the simulation (yr)
    time = np.arange(0, tmax + dt, dt)  # the time axis for the model
    dicts = add_emissions(copy_dicts(init_dicts), time, 600, 800, 8)

    vars_to_plot = ['DIC', 'TA', 'pCO2', 'particle_sinking_time', 'f_CaCO3']

    time_array, finished_dicts = original_model(copy_dicts(dicts), 2000, 0.5)
    final_lolat, final_hilat, final_deep, final_atmos = finished_dicts
    fig, axs = plot.boxes(time_array, vars_to_plot, final_lolat, final_hilat, final_deep, final_atmos,
                          label='Original Model')

    time_array, finished_dicts = acidification_model(copy_dicts(dicts), 2000, 0.5)
    final_lolat, final_hilat, final_deep, final_atmos = finished_dicts
    plot.boxes(time_array, vars_to_plot, final_lolat, final_hilat, final_deep, final_atmos,
               axs=axs, label='Acidification Model', ls='dotted')

    time_array, finished_dicts = ballasting_model(copy_dicts(dicts), 2000, 0.5)
    final_lolat, final_hilat, final_deep, final_atmos = finished_dicts
    plot.boxes(time_array, vars_to_plot, final_lolat, final_hilat, final_deep, final_atmos,
               axs=axs, label='Ballasting Feedback Model', ls='dashed')

    plt.suptitle('Plot of model variables for different models')
    fig.tight_layout()
    add_fancy_labels(axs, vars_to_plot)
    axs[-1].set_xlabel('Time (Years)')
    for ax in axs:
        maximum = ax.get_ylim()[1]
        minimum = ax.get_ylim()[0]
        ax.set_ylim(minimum, maximum)
        ax.fill_between((600, 800), (minimum, minimum), (maximum, maximum), color='lightgray')

    # plot the graphs
    plt.savefig('pics/LongQuestionPlot1.png', dpi=600)
    plt.show()


def plot_against_emissions(gauss=False):

    # this line of code runs initialise_dicts which sets up the initial dictionaries
    init_dicts = init_dicts_16()
    tmax = 2000  # how many years to simulate (yr)
    dt = 0.5  # the time step of the simulation (yr)
    time = np.arange(0, tmax + dt, dt)  # the time axis for the model
    dicts = add_emissions(copy_dicts(init_dicts), time, 600, 800, 8, gaussian=gauss)

    vars_to_plot = ['DIC', 'TA', 'pCO2', 'particle_sinking_time', 'f_CaCO3']
    # PLOT THE THREE MODELS INDIVIDUALLY COMPARED TO NO EMISSIONS

    models = [ballasting_model, acidification_model, original_model]
    names = ['Ballasting Feedback', 'Acidification', 'Original']
    for i, model in enumerate(models):

        # DECIDE WHICH VARS TO PLOT
        vars_to_plot = ['DIC', 'pCO2'] if names[i] == 'Original' else vars_to_plot
        vars_to_plot = ['DIC', 'TA', 'pCO2', 'f_CaCO3'] if names[i] == 'Acidification' else vars_to_plot
        gtype = ""
        if gauss:
            vars_to_plot.append('GtC_emissions')
            gtype = "(Gaussian)"

        # PLOT WITH EMISSIONS
        time_array, finished_dicts = model(copy_dicts(dicts), 2000, 0.5)
        final_lolat, final_hilat, final_deep, final_atmos = finished_dicts
        fig, axs = plot.boxes(time_array, vars_to_plot, final_lolat, final_hilat, final_deep, final_atmos,
                              label='Emissions', ls='solid')

        # PLOT WITHOUT EMISSIONS
        time_array, finished_dicts = model(copy_dicts(init_dicts), 2000, 0.5)
        final_lolat, final_hilat, final_deep, final_atmos = finished_dicts
        plot.boxes(time_array, vars_to_plot, final_lolat, final_hilat, final_deep, final_atmos,
                   axs=axs, label='No Emissions', ls='dotted')

        # FORMAT NICELY
        plt.suptitle(f'Plot of model variables for {names[i]} Model')
        fig.tight_layout()
        add_fancy_labels(axs, vars_to_plot)
        axs[-1].set_xlabel('Time (Years)')

        # ADD GREY EMISSION BARS
        if not gauss:
            for ax in axs:
                maximum = ax.get_ylim()[1]
                minimum = ax.get_ylim()[0]
                ax.set_ylim(minimum, maximum)
                ax.fill_between((600, 800), (minimum, minimum), (maximum, maximum), color='lightgray')

        plt.savefig(f'pics/LongQuestionPlot{names[i].replace(" ", "")}Model{gtype}.png', dpi=600)
        plt.show()


def plot_export():
    """
    this runs the code
    """

    # this line of code runs initialise_dicts which sets up the initial dictionaries
    init_dicts = init_dicts_16()
    tmax = 2000  # how many years to simulate (yr)
    dt = 0.5  # the time step of the simulation (yr)
    time = np.arange(0, tmax + dt, dt)  # the time axis for the model
    dicts = add_emissions(copy_dicts(init_dicts), time, 600, 800, 8)

    vars_to_plot = ['DIC_export', 'TA_export', 'f_CaCO3', 'exp']

    time_array, finished_dicts = ballasting_model(copy_dicts(dicts), 2000, 0.5)
    final_lolat, final_hilat, final_deep, final_atmos = finished_dicts
    fig, axs = plot.boxes(time_array, vars_to_plot, final_lolat, final_hilat,
                          ls='dashed', label='Ballasting Model')

    time_array, finished_dicts = acidification_model(copy_dicts(dicts), 2000, 0.5)
    final_lolat, final_hilat, final_deep, final_atmos = finished_dicts
    fig, axs = plot.boxes(time_array, vars_to_plot, final_lolat, final_hilat,
                          axs=axs, label='Acidification Model', ls='dotted')

    fig.tight_layout()
    add_fancy_labels(axs, vars_to_plot)
    axs[-1].set_xlabel('Time (Years)')

    # ADD GREY EMISSION BARS
    for ax in axs:
        maximum = ax.get_ylim()[1]
        minimum = ax.get_ylim()[0]
        ax.set_ylim(minimum, maximum)
        ax.fill_between((600, 800), (minimum, minimum), (maximum, maximum), color='lightgray')

    plt.savefig(f'pics/LongQuestionPlot2.png', dpi=600)
    plt.show()


if __name__ == '__main__':
    plot_export()
    # plot_all_three()
    # plot_against_emissions(gauss=False)
    # plot_against_emissions(gauss=True)
