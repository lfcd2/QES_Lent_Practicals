import matplotlib.pyplot as plt
import numpy as np

rho_org = 1200
rho_CaCO3 = 2700

lolat = {
    'name': 'Low Latitude',
    'rho_p':  1568.14,
    'depth': 100.,
    'k_ballast': 0.1609,
    'c': 'black'
}

hilat = {
    'name': 'High Latitude',
    'rho_p':  1416.27,
    'depth': 200.,
    'k_ballast': 0.0805,
    'c': 'red'
}

fCaCO3 = np.linspace(0, 1, 10000)
for box in [lolat, hilat]:
    rho_p = (rho_org + fCaCO3 * 10 / 3 * rho_org) / (1 + fCaCO3 * 10 / 3 * rho_org / rho_CaCO3)
    v = 10 * (rho_p - 1000) / (box['rho_p'] - 1000)
    box['particle_sinking_time'] = box['depth'] / v
    organic_export_efficiency = np.exp(- box['k_ballast'] * box['particle_sinking_time'])
    CaCO3_export_efficiency = fCaCO3 * organic_export_efficiency
    print(organic_export_efficiency, rho_p)
    plt.plot(fCaCO3, organic_export_efficiency,
             label=f'{box["name"]} organic', ls='solid', c=box['c'])
    plt.plot(fCaCO3, CaCO3_export_efficiency,
             label=f'{box["name"]} CaCO3', ls='dotted', c=box['c'])
plt.legend() # title='Export Efficiency')
plt.xlabel(r'$f_{CaCO_{3}}$')
plt.ylabel('Export Efficiency')
# plt.title(r'Plot of export efficiencies against $f_{CaCO_{3}}$')
plt.xlim(0, 1)
plt.tight_layout()
plt.savefig('pics/Efficiency and fCaCO3 plot.png', dpi=600)
plt.show()
