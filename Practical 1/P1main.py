import numpy as np
from matplotlib import pyplot as plt


def plot1():
    D = np.logspace(-2, 2, num=1000)

    IL = 300
    IH = 170
    B = 2.71
    DT = (IL - IH)/(4*D + B)
    Tav = 282.5
    TL = Tav + DT/2
    TH = Tav - DT/2
    F = 2*D*(TL - TH)

    Entropy = F/TH - F/TL

    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(10, 10))

    axs[0].plot(D,TH-273)
    axs[0].plot(D,TL-273)
    axs[1].plot(D,Entropy)
    axs[0].axvline(D[np.argmax(Entropy)],c='k',ls='--')

    for ax in axs:
        ax.set_xscale('log')
        ax.set_xlabel('Meridional heat diffusion coefficient')

    axs[0].set_ylabel('Temperature (°C)')
    axs[1].set_ylabel('Entropy (W$m^{-2}$$K^{-1}$)')

    axs[0].fill_between(D,27,15,color='C1', alpha=0.2)
    axs[0].fill_between(D,-5,5,color='C0', alpha=0.2)

    axs[1].text(s=f'D at MEP={D[np.argmax(Entropy)]:.2f}',x=0.75,y=0.9,transform=axs[1].transAxes)
    plt.show()


def plot2():
    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(10, 15), sharex=True)

    no_steps = 365
    f = np.zeros(no_steps)
    TL = np.zeros(no_steps)
    TH = np.zeros(no_steps)
    EL = np.zeros(no_steps)
    EH = np.zeros(no_steps)
    TL[0] = 298
    TH[0] = 280
    D = 1.35
    A = 203.3
    B = 2.71
    Ca = 1.02e7
    EL[0] = A + B * (TL[0] - 273)
    EH[0] = A + B * (TH[0] - 273)
    IL = 300
    IH = 170
    max_t = 365
    steps = np.linspace(0, max_t, no_steps)
    stepsize = max_t / no_steps

    for t in range(0, len(steps) - 1):

        EL[t+1] = A + B * (TL[t] - 273)
        EH[t+1] = A + B * (TH[t] - 273)
        TL[t+1] = TL[t] + (IL - EL[t] - 2 * D * (TL[t] - TH[t])) * (86400/Ca)
        TH[t+1] = TH[t] + (IH - EH[t] + 2 * D * (TL[t] - TH[t])) * (86400/Ca)
        f[t] = 2 * D * (TL[t] - TH[t])


    axs[0].plot(steps, TL - 273, 'b--', )
    axs[0].plot(steps, TH - 273, 'b-.')

    axs[1].plot(steps, EL, 'b--')
    axs[1].plot(steps, EH, 'b-.')

    axs[2].plot(steps[:-1], f[:-1])

    axs[1].plot([0, max_t], [IL, IL], 'k--')
    axs[1].plot([0, max_t], [IH, IH], 'k-.')
    axs[0].set_ylabel('Temperature ( ̊C)')
    axs[1].set_ylabel('Heat Flux (Wm$^{-2}$)')
    axs[2].set_ylabel('Latitudinal Heat Flow (Wm$^{-2}$)')
    axs[2].set_xlabel('Time (Days)')
    plt.show()


if __name__ == "__main__":
    plot1()
