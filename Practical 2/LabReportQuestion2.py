import numpy as np
import math
from matplotlib import pyplot as plt


def main():
    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(10, 15), sharex=True)
    for i, tau in enumerate([2, 10, 50, 100, 200]):
        # Set various model parameters
        alpha = 2e-4 # thermal expansion coefficient
        beta = 7e-4 # haline contraction coefficient
        k = 8.3e17 # proportionality coefficient to convert change in T and S to volume flux in units of m³/year
        # tau = 2 # Timescale for atmosphere/ocean temperature relaxation in years

        # Set the areas and volumes of the boxes
        SA_ocean = 3.58e14 # Surface area of the ocean
        AL = SA_ocean * 0.85 # Low latitude area units m²
        AH = SA_ocean * 0.15 # High latitude area units m²
        D = 3000 # depth of each box in meters
        VL = AL * D # volume of the low latitude box
        VH = AH * D # volume of the high latitude box
        Fw = 0.25 # high latitude evaporation minus precipitation in units of m/year
        Sref = 35 # reference salinity in units of g/kg (for calculation of E-P)
        E = AH * Fw * Sref # Evaporation minus precipitation (E-P)

        # Timestepping  parameters
        duration = 1000 # duration of the simulation in years
        steps = 1000 # Number of timesteps
        stepsize = duration/steps # timestep size
        times = np.linspace(0, duration, steps)

        # Create empty arrays for the low and high latitude T+S and volume flux
        TL = np.empty_like(times) # Initialize an array for the low latitude temperature
        SL = np.empty_like(times) # Initialize an array for the low latitude salinity
        TH = np.empty_like(times) # Initialize an array for the high latitude temperature
        SH = np.empty_like(times) # Initialize an array for the high latitude salinity
        Q = np.empty_like(times) # Initialize an array for the volume flux

        # Set initial conditions
        TL[0] = 298. # 298 low latitude temperature in degrees K
        SL[0] = 37. # low latitude salinity in PSS
        TH[0] = 273. # high latitude temperature in degrees K
        SH[0] = 33. # 33 high latitude salinity in PSS
        TatL = 20. + 273. # low latitude atmospheric temperature
        TatH = 0. + 273. # high latitude atmospheric temperature


        for t in range(0, steps-1):
            dT = TL[t] - TH[t]
            dS = SL[t] - SH[t]

            q = k*(alpha*dT - beta*dS)
            if t == 0:
                Q[t] = q
            Q[t+1] = q

            #  Part 2
            if 100 <= t <= 200:
                SH[t+1] = SH[t] - 0.000 + stepsize*(q*dS - E) / VH
            else:
                SH[t + 1] = SH[t] + stepsize * (q * dS - E) / VH

            TL[t+1] = TL[t] + stepsize*(-q*dT - (VL/tau)*(TL[t]-TatL)) / VL
            TH[t+1] = TH[t] + stepsize*(q*dT - (VH/tau)*(TH[t]-TatH)) / VH
            SL[t+1] = SL[t] + stepsize*(-q*dS + E) / VL


            # <Enter equations for calculating delT, delS, Q, and, if Q is greater than zero, TL, TH, SL and SH for each time step based on the values of the previous timestep>

        # calculate the circulation at the final time:
        # <Enter equations for calculating delT, delS, for the final timestep>

        # Create a figure to plot changes in temperature, salinity and volume flux

        if i == 0:
            line = axs[0].plot(times,TL-273,'--',label=r'Low Latitude')
            c = line[-1].get_color()
            axs[0].plot(times,TH-273,'-.',label=r'High Latitude', color=c)
        else:
            line = axs[0].plot(times, TL - 273, '--')
            c = line[-1].get_color()
            axs[0].plot(times, TH - 273, '-.', color=c)

        axs[1].plot(times,SL,'--',label=r'$S_L$', color=c)
        axs[1].plot(times,SH,'-.',label=r'$S_H$', color=c)

        axs[2].plot(times,Q, color=c, label=tau)
        #print(Q[-1], TL[-1], TH[-1], SH[-1], SL[-1])
        # Decorate figure
        axs[0].legend()
        axs[2].legend(title='Temperature\nRelaxation')
        axs[0].set_ylabel('Temperature ( ̊C)')
        axs[1].set_ylabel('Salinity')
        axs[2].set_ylabel('Latitudinal flow strength ($m^{3}yr^{-1}$)')
        axs[2].set_xlabel('Time (years)')
        # axs[0].set_ylim(0, 1)



if __name__ == "__main__":
    main()
    plt.show()