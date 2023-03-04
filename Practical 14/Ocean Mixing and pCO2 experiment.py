from P4main import initialise_dicts, ocean_model_p4
from OceanTools.tools import plot
from lfcd2OceanTools.lfcd2Tools import copy_dicts
import matplotlib.pyplot as plt


def experiment_ocean_mixing_and_pco2():
    # run the first model and plot it
    dicts = initialise_dicts()
    time_array, finished_dicts = ocean_model_p4(dicts, 3000, 0.5)
    lolat, hilat, deep, atmos = finished_dicts
    fig, axs = plot.boxes(time_array, 'pCO2', atmos)

    # make a copy and edit the dicts
    new_dicts = copy_dicts(dicts)
    for d in new_dicts:
        if d['name'] in ['hilat', 'lolat']:
            d['tau_M'] /= 2

    # run the updated model and plot it
    time_array, finished_dicts = ocean_model_p4(new_dicts, 3000, 0.5)
    lolat, hilat, deep, atmos = finished_dicts
    fig, axs = plot.boxes(time_array, 'pCO2', atmos, axs=axs, ls='--')

    # add legend (oscar's function seems broken for some reason?)
    lines = axs[0].get_lines()
    lines[1].set_label('Atmos')
    lines[2].set_label('Atmos: 2x Mixing')
    plt.legend()

    # show plot
    plt.show()


if __name__ == '__main__':
    experiment_ocean_mixing_and_pco2()
