from P5main import initialise_dicts_15, ocean_model_p5
from lfcd2OceanTools.lfcd2Tools import copy_dicts, modified_boxes
import numpy as np
from OceanTools.tools import plot
import matplotlib.pyplot as plt


def run():

    dicts = list(initialise_dicts_15())

    # create a new time axis for the model containing 3000 years with a 0.5 year time step
    tmax = 1000  # how many years to simulate (yr)
    dt = 0.5  # the time step of the simulation (yr)
    time = np.arange(0, tmax + dt, dt)  # the time axis for the model

    emit_atmos = dicts[-1].copy()  # create a copy of the original atmosphere input dictionary
    emit_atmos['GtC_emissions'] = np.zeros(time.shape)  # creat an array to hold the emission scenario
    emit_atmos['GtC_emissions'][(time > 200) & (time <= 400)] = 8.0  # set e to 8 GtC per year between 500-700
    dicts[-1] = emit_atmos

    acidification_dicts = copy_dicts(dicts)
    acidification_dicts[0]['f_CaCO3'] /= 2
    acidification_dicts[1]['f_CaCO3'] /= 2

    ballasting_dicts = copy_dicts(acidification_dicts)
    ballasting_dicts[0]['tau_PO4'] *= 2
    ballasting_dicts[1]['tau_PO4'] *= 2

    for i, a in enumerate([('solid', dicts), ('dotted', acidification_dicts), ('dashed', ballasting_dicts)]):
        ls, dicts2 = a
        time_array, finished_dicts = ocean_model_p5(dicts2, tmax, dt)
        final_lolat, final_hilat, final_deep, final_atmos = finished_dicts
        if i == 0:
            fig, axs = modified_boxes(time_array, ['pCO2'],
                                  final_lolat, final_hilat, final_deep, final_atmos, ls=ls, height=5)
        else:
            fig, axs = plot.boxes(time_array, ['pCO2'],
                                  final_lolat, final_hilat, final_deep, final_atmos, axs=axs, ls=ls)

    ax = axs[0]
    plt.plot([-1, -2], [0, 0], c='black', ls='solid', label='Original Model')
    plt.plot([-1, -2], [0, 0], c='black', ls='dotted', label='Ocean Acidification')
    plt.plot([-1, -2], [0, 0], c='black', ls='dashed', label='Ballasting Feedback')

    hand, leg = ax.get_legend_handles_labels()
    ax.legend(hand, leg)

    ax.set_ylim(0, 1200)
    ax.set_xlim(0, tmax)
    plt.show()


if __name__ == '__main__':
    run()
