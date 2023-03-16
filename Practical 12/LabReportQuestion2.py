"""
Some comments were added by Louis De Neve for the sake of poor QES students who are stuck
Finished, almost completely PEP-8 compliant code for most of the lab report is available on my GitHub:
https://github.com/lfcd2/QES_Lent_Practicals
"""

import numpy as np
from matplotlib import pyplot as plt


def main():
    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(8, 10), sharex=True)
    for i, tau in enumerate([2, 10, 25, 50, 100, 150, 200]):
        # Set various model parameters
        alpha = 2e-4  # thermal expansion coefficient
        beta = 7e-4  # haline contraction coefficient
        k = 8.3e17  # proportionality coefficient to convert change in T and S to volume flux in units of m³/year
        # tau = 2  # Timescale for atmosphere/ocean temperature relaxation in years

        # Set the areas and volumes of the boxes
        SA_ocean = 3.58e14  # Surface area of the ocean
        AL = SA_ocean * 0.85  # Low latitude area units m²
        AH = SA_ocean * 0.15  # High latitude area units m²
        D = 3000  # depth of each box in meters
        VL = AL * D  # volume of the low latitude box
        VH = AH * D  # volume of the high latitude box
        Fw = 0.25  # high latitude evaporation minus precipitation in units of m/year
        Sref = 35  # reference salinity in units of g/kg (for calculation of E-P)
        E = AH * Fw * Sref  # Evaporation minus precipitation (E-P)

        # Timestepping  parameters
        duration = 1000  # duration of the simulation in years
        steps = 10000  # Number of timesteps
        stepsize = duration / steps  # timestep size
        times = np.linspace(0, duration, steps)

        # Create empty arrays for the low and high latitude T+S and volume flux
        TL = np.empty_like(times)  # Initialize an array for the low latitude temperature
        SL = np.empty_like(times)  # Initialize an array for the low latitude salinity
        TH = np.empty_like(times)  # Initialize an array for the high latitude temperature
        SH = np.empty_like(times)  # Initialize an array for the high latitude salinity
        Q = np.empty_like(times)  # Initialize an array for the volume flux

        # Set initial conditions
        TL[0] = 298.  # low latitude temperature in degrees K
        SL[0] = 37.  # low latitude salinity in PSS
        TH[0] = 273.  # high latitude temperature in degrees K
        SH[0] = 33.  # high latitude salinity in PSS
        TatL = 20. + 273.  # low latitude atmospheric temperature
        TatH = 0. + 273.  # high latitude atmospheric temperature

        for t in range(0, steps-1):

            # calculates the deltas and q
            dT = TL[t] - TH[t]
            dS = SL[t] - SH[t]
            q = k * (alpha * dT - beta * dS)

            # assigns the current q to the array
            Q[t] = q

            # uses the differential equations given to calculate the next TL, TH, SL and SH
            TL[t + 1] = TL[t] + stepsize * (-q * dT - (VL / tau) * (TL[t] - TatL)) / VL
            TH[t + 1] = TH[t] + stepsize * (q * dT - (VH / tau) * (TH[t] - TatH)) / VH
            SL[t + 1] = SL[t] + stepsize * (-q * dS + E) / VL
            SH[t + 1] = SH[t] + stepsize * (q * dS - E) / VH

        # this bit of code calculates the final value of Q to prevent it from being = 0
        dT = TL[-1] - TH[-1]
        dS = SL[-1] - SH[-1]
        Q[-1] = k * (alpha * dT - beta * dS)

        # This plots all the suitable lines, with matching colours. It only adds a label if it is the first iteration
        if i == 0:
            line = axs[0].plot(times, TL-273, '--', label=r'Low Latitude')  # plots the first line
            c = line[-1].get_color()  # copies the colour from the first line
            axs[0].plot(times, TH-273, '-.', label=r'High Latitude', color=c)  # plots the matching line
        else:
            line = axs[0].plot(times, TL - 273, '--')  # plots the first line
            c = line[-1].get_color()  # copies the colour from the first line
            axs[0].plot(times, TH - 273, '-.', color=c)  # plots the matching line

        # plots the salinity line
        axs[1].plot(times, SL, '--', label=r'$S_L$', color=c)
        axs[1].plot(times, SH, '-.', label=r'$S_H$', color=c)

        # plots Q
        axs[2].plot(times, Q, color=c, label=tau)

    # Decorate figure
    axs[0].legend(loc='upper right')
    axs[2].legend(title='Temperature\nRelaxation', loc='upper right')
    axs[0].set_ylabel('Temperature ( ̊C)')
    axs[1].set_ylabel('Salinity (PSS)')
    axs[2].set_ylabel('Latitudinal flow strength ($m^{3}yr^{-1}$)')
    axs[2].set_xlabel('Time (years)')

    for ax in axs:
        ax.set_xlim(0, 1000)
    # show figure
    fig.tight_layout()
    plt.savefig('Lab_Report_Q2_Output.png', dpi=600)
    plt.show()


if __name__ == "__main__":
    main()
