import numpy as np
from matplotlib import pyplot as plt


# this function does the whole question
def plot_graph():

    # generates the axes and figure objects
    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(10, 12), sharex=True)
    # gives it a nice tite :)
    plt.suptitle('North and South hemisphere temperatures over time dependent on'
                 '\ninitial temperature and meridional heat diffusion coefficient', y=0.99)

    # declares a bunch of variables outside the loop for huge optimisation (edited from oscars code)
    no_steps = 365
    IL = 300
    IH = 170
    max_t = 365
    A = 202
    B = 2.71
    Ca = 1.02e7
    steps = np.linspace(0, max_t, no_steps)  # this step calculation is the same as in main part of practical
    stepsize = max_t / no_steps

    # okay now its actually almost getting complex
    # This line just means we're going to plot the same thing for 5 evenly spaced values of D
    for i, D in enumerate(np.linspace(0.35, 1.35, 5)):
        print(D)  # unnecessary print here but once again, we move

        # now, for each value of D, we will plit the same line three times for different starting values of temperature
        # the for x, y in [[a, b], [c, d]] just means repeat it, but have x = a, y = b the first time, and then
        # repeat for x = c, y = d.
        for tl, th in [[305, 283], [300, 278], [295, 273]]:

            # make some empty arrays of a suitable size
            f = np.zeros(no_steps)
            TL = np.zeros(no_steps)
            TH = np.zeros(no_steps)
            EL = np.zeros(no_steps)
            EH = np.zeros(no_steps)

            # assign starting values for T and E
            TL[0] = tl
            TH[0] = th
            EL[0] = A + B * (TL[0] - 273)
            EH[0] = A + B * (TH[0] - 273)

            # this is now the same euler method as prior, w/ the same equations (although written a bit differently, mb)
            for t in range(0, len(steps) - 1):
                EL[t + 1] = A + B * (TL[t] - 273)
                EH[t + 1] = A + B * (TH[t] - 273)
                TL[t + 1] = TL[t] + (IL - EL[t] - 2 * D * (TL[t] - TH[t])) * (86400 / Ca) * stepsize
                TH[t + 1] = TH[t] + (IH - EH[t] + 2 * D * (TL[t] - TH[t])) * (86400 / Ca) * stepsize
                f[t] = 2 * D * (TL[t] - TH[t])

            # plots temps
            axs[0].plot(steps, TL - 273, '--', c='C' + str(i))
            axs[0].plot(steps, TH - 273, '-.', c='C' + str(i))

            # plots Es
            axs[1].plot(steps, EL, '--', c='C' + str(i))
            axs[1].plot(steps, EH, '-.', c='C' + str(i))

            # plots the final latitudinal heat flow, but if it's the first iteration, also add a label
            if tl == 305:
                axs[2].plot(steps[:-1], f[:-1], c='C' + str(i), label=D)
            else:
                axs[2].plot(steps[:-1], f[:-1], c='C' + str(i))

    # we're now outside the loop so no longer are repeating everything

    # sets the x axes to go no further than max_t
    for ax in axs:
        ax.set_xlim(0, max_t)

    # adds reference bars, with labels, and then adds a legend
    axs[1].plot([0, max_t], [IL, IL], 'k--', label='Incoming Radiation (Low Latitude)')
    axs[1].plot([0, max_t], [IH, IH], 'k-.', label='Incoming Radiation (High Latitude)')
    axs[1].legend()

    # adds more text and labels
    axs[0].text(10, 15.5, 'Top sets = low latitudes, bottom sets = high latitudes, '
                          'Shaded areas represent range of normal temperatures', c='grey')
    axs[0].set_ylabel('Temperature ( ̊C)')
    axs[1].set_ylabel('Heat Flux (Wm$^{-2}$)')
    axs[2].set_ylabel('Latitudinal Heat Flow (Wm$^{-2}$)')
    axs[2].set_xlabel('Time (Days)')

    # adds the shaded areas, but makes them grey for clarity
    axs[0].fill_between(steps, 27, 15, color='lightgray')
    axs[0].fill_between(steps, -5, 5, color='lightgray')

    # adds the bottom legend, but with more fancy positioning etc
    axs[2].legend(loc='upper center', bbox_to_anchor=(0.761, 1.035), title='Value of D', ncol=3, fancybox=True)

    # makes it pretty :)
    fig.tight_layout()

    # saves the graph and displays it
    plt.savefig('Lab_Report_Q1_Output.png', dpi=600)
    plt.show()


# runs the code y'know
if __name__ == "__main__":
    plot_graph()
