"""
Many of the comments were added by Louis De Neve for the sake of poor QES students who are stuck
Finished, almost completely PEP-8 compliant code for most of the lab report is available on my GitHub:
https://github.com/lfcd2/QES_Lent_Practicals
"""

from P5main import initialise_dicts_15, ocean_model_p15
from lfcd2OceanTools.lfcd2Tools import copy_dicts, modified_boxes
import numpy as np
from OceanTools.tools import plot
import matplotlib.pyplot as plt


def run():
    # make a set of dictionaries, same way as in P5 main
    dicts = list(initialise_dicts_15())

    # create a new time axis for the model containing 3000 years with a 0.5 year time step
    tmax = 3000  # how many years to simulate (yr)
    dt = 0.5  # the time step of the simulation (yr)
    time = np.arange(0, tmax + dt, dt)  # the time axis for the model

    emit_atmos = dicts[-1].copy()  # create a copy of the original atmosphere input dictionary
    emit_atmos['GtC_emissions'] = np.zeros(time.shape)  # creat an array to hold the emission scenario
    emit_atmos['GtC_emissions'][(time > 200) & (time <= 400)] = 8.0  # set e to 8 GtC per year between 500-700
    dicts[-1] = emit_atmos  # insert the modified dict back into dicts

    acidification_dicts = copy_dicts(dicts)  # makes a copy of dicts and changes f_CaCO3
    acidification_dicts[0]['f_CaCO3'] /= 2
    acidification_dicts[1]['f_CaCO3'] /= 2

    ballasting_dicts = copy_dicts(acidification_dicts)  # makes a copy of the new dicts and increases tau_PO4
    ballasting_dicts[0]['tau_PO4'] *= 2
    ballasting_dicts[1]['tau_PO4'] *= 2

    # this iterates over the three model dicts
    for i, a in enumerate([('solid', 'Original Model', dicts),
                           ('dotted', 'Ocean Acidification', acidification_dicts),
                           ('dashed', 'Ballasting Feedback', ballasting_dicts)]):
        ls, name, dicts2 = a

        # for every set of intial conditions, this runs the model and unpacks the dicts
        time_array, finished_dicts = ocean_model_p15(dicts2, tmax, dt)
        final_lolat, final_hilat, final_deep, final_atmos = finished_dicts

        # if this is the model, it will make a new fig, axs and plot the graph
        if i == 0:
            fig, axs = modified_boxes(time_array, ['pCO2'],
                                      final_lolat, final_hilat, final_deep, final_atmos, ls=ls, height=5)
        # if it's the second or third, this is the same, but it inputs the axes
        else:
            fig, axs = plot.boxes(time_array, ['pCO2'],
                                  final_lolat, final_hilat, final_deep, final_atmos, axs=axs, ls=ls)

        # plots some hidden, very small lines to add to the legend
        plt.plot([-1, -2], [0, 0], c='black', ls=ls, label=name)

        # prints the equilibrium values
        print('Final pCO2 Values:')
        for box in [final_lolat, final_hilat, final_atmos]:
            print(f'Model: {name}, box: {box["name"]} : {box["pCO2"][-1]} ')

    # makes a legend and adjusts limits
    ax = axs[0]
    ax.legend()
    ax.set_ylim(0, 1200)
    ax.set_xlim(0, tmax)

    # plots it
    plt.show()


if __name__ == '__main__':
    run()
