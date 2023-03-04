"""
Some comments were added by Louis De Neve for the sake of poor QES students who are stuck
Finished, almost completely PEP-8 compliant code for most of the lab report is available on my GitHub:
https://github.com/lfcd2/QES_Lent_Practicals
"""

import numpy as np
from matplotlib import pyplot as plt


def main():
    # Set various model parameters
    alpha = 2e-4  # thermal expansion coefficient
    beta = 7e-4  # haline contraction coefficient
    k = 8.3e17  # proportionality coefficient to convert change in T and S to volume flux in units of m³/year
    tau = 2  # Timescale for atmosphere/ocean temperature relaxation in years

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
    steps = 1000  # Number of timesteps
    stepsize = duration/steps  # timestep size
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

    # these are the initial conditions needed for part 2
    # TL[0] = 292.87
    # TH[0] = 273.73
    # SH[0] = 36.2706
    # SL[0] = 36.4228

    # repeats for the number of steps
    for t in range(0, steps-1):

        # calculates the deltas and q
        dT = TL[t] - TH[t]
        dS = SL[t] - SH[t]
        q = k*(alpha*dT - beta*dS)

        # assigns the current q to the array
        Q[t] = q

        # uses the differential equations given to calculate the next TL, TH, SL and SH
        TL[t + 1] = TL[t] + stepsize * (-q * dT - (VL / tau) * (TL[t] - TatL)) / VL
        TH[t + 1] = TH[t] + stepsize * (q * dT - (VH / tau) * (TH[t] - TatH)) / VH
        SL[t + 1] = SL[t] + stepsize * (-q * dS + E) / VL
        SH[t + 1] = SH[t] + stepsize * (q * dS - E) / VH

        # this is part 2 - it decreases some salinity each cycle
        if 100 <= t <= 200:
            SH[t+1] = SH[t+1] - 0.002

    # this bit of code calculates the final value of Q to prevent it from being = 0
    dT = TL[-1] - TH[-1]
    dS = SL[-1] - SH[-1]
    Q[-1] = k * (alpha * dT - beta * dS)

    # Create a figure to plot changes in temperature, salinity and volume flux
    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(8, 10), sharex=True)

    # plots the temps
    axs[0].plot(times, TL-273, '--', label=r'$T_L$')
    axs[0].plot(times, TH-273, '-.', label=r'$T_H$')

    # plots the salinities
    axs[1].plot(times, SL, '--', label=r'$S_L$')
    axs[1].plot(times, SH, '-.', label=r'$S_H$')

    # plots q
    axs[2].plot(times, Q)
    print(Q[-1], TL[-1], TH[-1], SH[-1], SL[-1])

    # Decorate figure
    axs[0].legend()
    axs[1].legend()
    axs[0].set_ylabel('Temperature ( ̊C)')
    axs[1].set_ylabel('Salinity')
    axs[2].set_ylabel('Latitudinal flow strength ($m^{3}yr^{-1}$)')
    axs[2].set_xlabel('Time (years)')

    # shows the figure
    plt.show()


if __name__ == "__main__":
    main()
