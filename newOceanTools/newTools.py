import matplotlib.pyplot as plt
import numpy as np


def modified_boxes(time, vars, *boxes, axs=None, label=None, height=0, **kwargs):
    """Plot a set of variables in a set of boxes.

    Parameters
    ----------
    time : array-like
        An array containing the time axis of the model.
    vars : list of str
        A list containing the names of the variables to plot.
    *boxes : dicts
        The boxes that you want to plot. Each box is a dictionary.
    axs : list of Axes, optional
        A list of axes to plot the variables in. If not given, a new figure and axes are created.
    label : str, optional
        Some custom text to add to the legend.
    height : float, optional
        Height modifier to edit the height of the figure
    **kwargs
        Keyword arguments to pass to the pyplot.plot function.
    """
    if isinstance(vars, str):
        vars = [vars]

    if axs is None:
        fig, axs = plt.subplots(len(vars), 1, figsize=(8, (1.7 * len(vars)) + height), sharex=True, constrained_layout=True)
    else:
        fig = axs[0].figure

    if hasattr(axs, 'plot'):
        axs = [axs]

    cdict = {
        'deep': 'tab:grey',
        'hilat': 'tab:blue',
        'lolat': 'tab:red',
        'atmos': 'tab:orange',
    }

    for var, ax in zip(vars, axs):
        for box in boxes:
            if var in box:
                boxname = box['name']
                try:
                    ax.plot(time, box[var], color=cdict[boxname], **kwargs)
                except ValueError:
                    continue

        ax.set_ylabel(var)

    # box labels, if they're not already there
    current_labels = axs[-1].get_legend_handles_labels()[1]
    plot_orig = False
    for box in boxes:
        if box['name'] not in current_labels:
            axs[-1].plot([], [], color=cdict[box['name']], label=box['name'])
        else:
            plot_orig = True
    axs[-1].legend(fontsize=8)

    if (label is not None) and len(axs) > 1:
        current_labels = axs[-2].get_legend_handles_labels()[1]
        if plot_orig and 'original' not in current_labels:
            axs[-2].plot([], [], color=(.3, .3, .3), label='original')
        axs[-2].plot([], [], color=(.3, .3, .3), label=label, **kwargs)
        axs[-2].legend(fontsize=8)

    axs[-1].set_xlabel('time')
    axs[-1].set_xlim(0, max(time))

    return fig, axs


def copy_dicts(dicts):
    """Returns a copy of each dictionary in the list (or of an individual dictionary)

    Parameters
    ----------
    dicts : array-like or dict
        array of dictionaries to be copied

    Returns
    -------
    new_dicts : list
        list of copied dictionaries
    """

    if type(dicts) == dict:
        return dicts.copy()

    else:
        try:
            new_dicts = []
            for dictionary in dicts:
                new_dicts.append(dictionary.copy())

            return new_dicts

        except AttributeError:
            print('Dictionaries could not be copied, output was the original dicts')
            return dicts


def modify_single_dict(d, modifiers):
    """ Modifies some variables in a dictionary

    Parameters
    ----------
    d : dict
        Dictionary to be modified
    modifiers : Modifier or array-like
        Modifier or list of modifiers to be applied to the dictionaries

    Returns
    ----------
    d : dict
        return dictionary with modifers applied
    """

    if type(modifiers) == Modifier:
        d = modifiers.apply(d)
    else:
        for modifier in modifiers:
            d = modifier.apply(d)
    return d


def modify_dicts(dicts, modifiers):
    """ Modifies some variables in a list of dicts

        Parameters
        ----------
        dicts : dict or array-like
            Dictionary or list of dictionaries to be modified
        modifiers : Modifier or array-like
            Modifier or list of modifiers to be applied to the dictionaries

        Returns
        ----------
        d : dict or array-like
            returns modified dictionaries as individual dict or list of dicts
        """
    if type(dicts) == dict:
        new_dicts = modify_single_dict(dicts, modifiers)

    else:
        new_dicts = []
        for d in dicts:
            try:
                d = modify_single_dict(d, modifiers)
            except AttributeError:
                print(f'Could Not Modify Dictionary {d["name"]}, due to Attribute Error')
            new_dicts.append(d)

    return new_dicts


class Modifier:
    def __init__(self,
                 boxes,
                 variables,
                 coeff=1,
                 multiply=True
                 ):
        self.variables = variables
        self.boxes = boxes
        self.coeff = coeff
        self.multiply = multiply

    def apply(self, d):
        if d['name'] in self.boxes:
            for var in self.variables:
                try:
                    if self.multiply:
                        d[var] *= self.coeff
                    else:
                        d[var] += self.coeff
                except TypeError:
                    print('Incorrect type modified')
        return d


def add_emissions(dicts, time, start, stop, value, gaussian=False):
    for i, d in enumerate(dicts):
        if d['name'] == 'atmos':
            emit_atmos = d.copy()  # create a copy of the original atmosphere input dictionary
            if type(emit_atmos['GtC_emissions']) == float:
                emit_atmos['GtC_emissions'] = np.zeros(time.shape)  # creat an array to hold the emission scenario
            if gaussian:
                SD = stop - start
                center = stop + start
                A = (value / np.sqrt(2*np.pi)) * 2
                for j, a in enumerate(emit_atmos['GtC_emissions']):
                    emit_atmos['GtC_emissions'][j] = A*np.exp(-((j-center)**2) / (2 * SD ** 2))
            else:
                emit_atmos['GtC_emissions'][(time > start) & (time <= stop)] = value  # set emissions
            dicts[i] = emit_atmos  # insert the modified dict back into dicts
    return dicts


def add_fancy_labels(axs, vars_to_plot):
    if 'f_CaCO3' in vars_to_plot:
        axs[vars_to_plot.index('f_CaCO3')].set_ylabel(r'$f_{CaCO_{3}}$ (unitless)')
    if 'pCO2' in vars_to_plot:
        axs[vars_to_plot.index('pCO2')].set_ylabel(r'$pCO_{2} (ppm)$')
    if 'DIC' in vars_to_plot:
        axs[vars_to_plot.index('DIC')].set_ylabel(r'DIC (mol m$^{-3}$)')
    if 'TA' in vars_to_plot:
        axs[vars_to_plot.index('TA')].set_ylabel(r'TA (mol m$^{-3}$)')
    if 'particle_sinking_time' in vars_to_plot:
        axs[vars_to_plot.index('particle_sinking_time')].set_ylabel(r'Sinking Time (days)')
    if 'GtC_emissions' in vars_to_plot:
        axs[vars_to_plot.index('GtC_emissions')].set_ylabel(r'Emissions (GtC)')
    if 'DIC_export' in vars_to_plot:
        axs[vars_to_plot.index('DIC_export')].set_ylabel(r'DIC Export (mol m$^{-3}$ yr$^{-1}$)', fontsize=8)
    if 'TA_export' in vars_to_plot:
        axs[vars_to_plot.index('TA_export')].set_ylabel(r'TA Export (mol m$^{-3}$ yr$^{-1}$)', fontsize=8)
    return axs
