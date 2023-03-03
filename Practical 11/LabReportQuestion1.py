import numpy as np
from matplotlib import pyplot as plt


def plot_graph():
    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(10, 12), sharex=True)
    plt.suptitle('North and South hemisphere temperatures over time dependent on'
                 '\ninitial temperature and meridional heat diffusion coefficient', y=0.99)
    cols = ['r', 'm', 'g', 'c', 'b']

    for i, D in enumerate(np.linspace(0.35, 1.35, 5)):
        print(D)
        for tl, th in [[305, 283], [300, 278], [295, 273]]:
            no_steps = 365
            f = np.zeros(no_steps)
            TL = np.zeros(no_steps)
            TH = np.zeros(no_steps)
            EL = np.zeros(no_steps)
            EH = np.zeros(no_steps)
            TL[0] = tl
            TH[0] = th

            A = 202
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
                EL[t + 1] = A + B * (TL[t] - 273)
                EH[t + 1] = A + B * (TH[t] - 273)
                TL[t + 1] = TL[t] + (IL - EL[t] - 2 * D * (TL[t] - TH[t])) * (86400 / Ca)
                TH[t + 1] = TH[t] + (IH - EH[t] + 2 * D * (TL[t] - TH[t])) * (86400 / Ca)
                f[t] = 2 * D * (TL[t] - TH[t])

            axs[0].plot(steps, TL - 273, '--', c='C' + str(i))
            axs[0].plot(steps, TH - 273, '-.', c='C' + str(i))

            axs[1].plot(steps, EL, '--', c='C' + str(i))
            axs[1].plot(steps, EH, '-.', c='C' + str(i))

            if tl == 305:
                axs[2].plot(steps[:-1], f[:-1], c='C' + str(i), label=D)
            else:
                axs[2].plot(steps[:-1], f[:-1], c='C' + str(i))

    for ax in axs:
        ax.set_xlim(0, max_t)

    axs[1].plot([0, max_t], [IL, IL], 'k--', label='Incoming Radiation (Low Latitude)')
    axs[1].plot([0, max_t], [IH, IH], 'k-.', label='Incoming Radiation (High Latitude)')

    axs[1].legend()
    axs[0].text(10, 15.5, 'Top sets = low latitudes, bottom sets = high latitudes, '
                          'Shaded areas represent range of normal temperatures', c='grey')
    axs[0].set_ylabel('Temperature ( ̊C)')
    axs[1].set_ylabel('Heat Flux (Wm$^{-2}$)')
    axs[2].set_ylabel('Latitudinal Heat Flow (Wm$^{-2}$)')
    axs[2].set_xlabel('Time (Days)')
    axs[0].fill_between(steps, 27, 15, color='lightgray')
    axs[0].fill_between(steps, -5, 5, color='lightgray')
    axs[2].legend(loc='upper center', bbox_to_anchor=(0.761, 1.035), title='Value of D', ncol=3, fancybox=True)

    fig.tight_layout()
    plt.savefig('Lab_Report_Q1_Output.png', dpi=300)
    plt.show()


if __name__ == "__main__":
    plot_graph()
