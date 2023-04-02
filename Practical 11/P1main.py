import numpy as np
from matplotlib import pyplot as plt


# first part of the practical
def plot1():

    # Initiate variables that will be used throughout
    D = np.logspace(-2, 2, num=1000)
    IL = 300
    IH = 170
    B = 2.71
    DT = (IL - IH)/(B + 4 * D)
    Tav = 282.5
    TL = Tav + DT/2
    TH = Tav - DT/2
    F = 2 * D * (TL - TH)
    Entropy = F/TH - F/TL

    # make figure and axes object with suitable rows etc
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(8, 10))

    # plots the high and low latitude temperatures (converting to Celsius with -273)
    axs[0].plot(D, TH-273)
    axs[0].plot(D, TL-273)

    # plots the entropy of the system on the second plot
    axs[1].plot(D, Entropy)

    # plots a vertical dashed line at the highest peak in the entropy curve
    axs[0].axvline(D[np.argmax(Entropy)], c='k', ls='--')

    # makes the fancy log scale on the x-axis
    for ax in axs:
        ax.set_xscale('log')
        ax.set_xlabel('Meridional heat diffusion coefficient')

    # relabels the y axes
    axs[0].set_ylabel('Temperature (°C)')
    axs[1].set_ylabel('Entropy (W$m^{-2}$$K^{-1}$)')

    # adds some boxes on the first plot to represent average ocean temperatures
    axs[0].fill_between(D, 27, 15, color='C1', alpha=0.2)
    axs[0].fill_between(D, -5, 5, color='C0', alpha=0.2)

    # adds some text with the x coordinate of the dotted line plotted earlier
    axs[1].text(s=f'D at MEP={D[np.argmax(Entropy)]:.2f}', x=0.75, y=0.9, transform=axs[1].transAxes)
    plt.show()


# second part of the practical
def plot2():

    # makes figure and axes objects
    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(8, 10), sharex=True)

    # creates empty arrays of suitable size for f, TL, TH, EL, EH, and sets starting values of T
    no_steps = 365
    f = np.zeros(no_steps)
    TL = np.zeros(no_steps)
    TH = np.zeros(no_steps)
    EL = np.zeros(no_steps)
    EH = np.zeros(no_steps)
    TL[0] = 298
    TH[0] = 280

    # defines a bunch more variables
    D = 0.85
    A = 203.3
    B = 2.71
    Ca = 1.02e7

    # sets starting values of E
    EL[0] = A + B * (TL[0] - 273)
    EH[0] = A + B * (TH[0] - 273)

    # declares even more variables (why didn't oscar just do this at the start of the code IDK im kinda done)
    IL = 300
    IH = 170
    max_t = 365

    # calculates the size of each increment (stepsize) based off the number of steps and timescale
    steps = np.linspace(0, max_t, no_steps)
    stepsize = max_t / no_steps

    # iterates for each time step (again, this is a really stupid way of doing it, but we move)
    for t in range(0, len(steps) - 1):

        # these are the  solutions of the differential equations given where stepsize * 86400 is t in seconds
        EL[t+1] = A + B * (TL[t] - 273)
        EH[t+1] = A + B * (TH[t] - 273)
        TL[t+1] = TL[t] + (IL - EL[t] - 2 * D * (TL[t] - TH[t])) * (stepsize*86400/Ca)
        TH[t+1] = TH[t] + (IH - EH[t] + 2 * D * (TL[t] - TH[t])) * (stepsize*86400/Ca)
        f[t] = 2 * D * (TL[t] - TH[t])
        # i won't explain these too much as it's the whole objective of the practical

    # plots the temperature on the first axes
    axs[0].plot(steps, TL - 273, 'b--', )
    axs[0].plot(steps, TH - 273, 'b-.')

    # plots heat flux on the second axes
    axs[1].plot(steps, EL, 'b--')
    axs[1].plot(steps, EH, 'b-.')

    # plots latitudinal heat flow on final axes
    axs[2].plot(steps[:-1], f[:-1])

    # adds some reference lines to the plots
    axs[1].plot([0, max_t], [IL, IL], 'k--')
    axs[1].plot([0, max_t], [IH, IH], 'k-.')

    # labels
    axs[0].set_ylabel('Temperature ( ̊C)')
    axs[1].set_ylabel('Heat Flux (Wm$^{-2}$)')
    axs[2].set_ylabel('Latitudinal Heat Flow (Wm$^{-2}$)')
    axs[2].set_xlabel('Time (Days)')

    # shows the plot
    plt.show()


# runs the code y'know
if __name__ == "__main__":
    plot1()
    plot2()
