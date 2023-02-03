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



if __name__ == "__main__":
  plot1()