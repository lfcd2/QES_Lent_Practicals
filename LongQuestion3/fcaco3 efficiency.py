import matplotlib.pyplot as plt
import numpy as np

rho_org = 1100
rho_CaCO3 = 2700

lolat = {
    'name': 'lolat',
    'rho_p':  1568.14,
    'depth': 100.,
    'k_ballast': 0.1609
}

hilat = {
    'name': 'hilat',
    'rho_p':  1416.27,
    'depth': 200.,
    'k_ballast': 0.0805
}

fCaCO3 = np.linspace(0, 1, 10000)
for box in [lolat, hilat]:
    rho_p = (rho_org + fCaCO3 * 10 / 3 * rho_org) / (1 + fCaCO3 * 10 / 3 * rho_org / rho_CaCO3)
    v = 10 * (rho_p - 1000) / (box['rho_p'] - 1000)
    box['particle_sinking_time'] = box['depth'] / v
    exponential_factor = np.exp(- box['k_ballast'] * box['particle_sinking_time'])
    print(exponential_factor, rho_p)
    plt.plot(fCaCO3, exponential_factor)
plt.xlabel('f_CaCO3')
plt.ylabel('exponential')
plt.show()
