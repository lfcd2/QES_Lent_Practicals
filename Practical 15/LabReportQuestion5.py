"""
Many of the comments were added by Louis De Neve for the sake of poor QES students who are stuck
Finished, almost completely PEP-8 compliant code for most of the lab report is available on my GitHub:
https://github.com/lfcd2/QES_Lent_Practicals
"""

from P5main import initialise_dicts_15, ocean_model_p15
from lfcd2OceanTools.lfcd2Tools import copy_dicts, Modifier, modify_dicts, add_emissions
import numpy as np
from OceanTools.tools import plot
import matplotlib.pyplot as plt


def run():
    # make a set of dictionaries, same way as in P5 main
    dicts = list(initialise_dicts_15())

    # create a new time axis for the model containing 3000 years with a 0.5 year time step
    tmax = 2000  # how many years to simulate (yr)
    dt = 0.5  # the time step of the simulation (yr)
    time = np.arange(0, tmax + dt, dt)  # the time axis for the model

    dicts = add_emissions(dicts, time, 200, 400, 8)  # adds 8GtC emissions for 200-400

    m = Modifier(['hilat', 'lolat'], ['f_CaCO3'], 0.5)
    acidification_dicts = modify_dicts(copy_dicts(dicts), m)

    m = Modifier(['hilat', 'lolat'], ['tau_PO4'], 2)
    ballasting_dicts = modify_dicts((copy_dicts(acidification_dicts)), m)
    # makes a copy of the new dicts and increases tau_PO4

    # this iterates over the three model dicts
    for i, a in enumerate([('solid', 'Original Model', dicts),
                           ('dotted', 'Ocean Acidification', acidification_dicts),
                           ('dashed', 'Ballasting Feedback', ballasting_dicts)]):
        ls, name, dicts2 = a

        # for every set of initial conditions, this runs the model and unpacks the dicts
        time_array, finished_dicts = ocean_model_p15(dicts2, tmax, dt)
        final_lolat, final_hilat, final_deep, final_atmos = finished_dicts

        # if this is the model, it will make a new fig, axs and plot the graph
        if i == 0:
            fig, axs = plot.boxes(time_array, ['pCO2', 'DIC', 'TA'],
                                  final_lolat, final_hilat, final_deep, final_atmos, ls=ls, height=2.5)
        # if it's the second or third, this is the same, but it inputs the axes
        else:
            fig, axs = plot.boxes(time_array, ['pCO2', 'DIC', 'TA'],
                                  final_lolat, final_hilat, final_deep, final_atmos, axs=axs, ls=ls, label=name)

        # prints the equilibrium values
        print('Final pCO2 Values:')
        for box in [final_atmos]:
            print(f'Model: {name}, box: {box["name"]} : {box["pCO2"][-1]} ')

    # adjusts limits
    ax = axs[0]
    ax.set_ylim(200, 1150)
    ax.set_xlim(0, tmax)
    ax.set_ylabel(r'pCO$_2$ (ppm)')
    axs[2].set_ylabel(r'TA (mol m$^{-3}$)')
    axs[1].set_ylabel(r'DIC (mol m$^{-3}$)')
    axs[2].set_xlabel('Time (Years)')
    axs[1].legend(loc=(0.7, 0.6))
    axs[2].legend(loc='lower right')
    # plt.suptitle('Plot of model variables against time')
    fig.tight_layout()

    for ax in axs:
        maximum = ax.get_ylim()[1]
        minimum = ax.get_ylim()[0]
        ax.set_ylim(minimum, maximum)
        ax.fill_between((200, 400), (minimum, minimum), (maximum, maximum), color='lightgray')

    # plots it
    plt.savefig('Lab_Report_Q5_Output.png', dpi=600)
    plt.show()


if __name__ == '__main__':
    run()
