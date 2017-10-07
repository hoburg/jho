import numpy as np
import matplotlib.pyplot as plt

Nmax = 3.5
N = np.linspace(0, Nmax, 500)
mtow = 713.75 # N
S = 2.0903 # m^2b
rho = 0.77070838 # kg/m^3
rho = 1.225 # kg/m^3
CLmax = 1.39

vstall = (2*mtow*N/S/rho/CLmax)**0.5
vmax = 44.01 # m/s

fig, ax = plt.subplots()
ax.plot(vstall, N, 'k')
ax.plot([max(vstall), vmax], [Nmax, Nmax], 'k')
ax.plot([vmax, vmax], [-2, 3.5], 'k')
N = np.linspace(0, 2, 500)
vstall = (2*mtow*N/S/rho/CLmax)**0.5
ax.plot([max(vstall), vmax], [-2, -2], 'k')
ax.plot(vstall, -N, 'k')
ax.set_xlim([0, 55])
ax.set_ylim([-2.5, 5])
ax.text(36, 3.9, 'max g-load',)
ax.text(47, 0.5, 'max speed',)
ax.set_xlabel("airspeed [m/s]")
ax.set_ylabel("load factor")
ax.grid()
fig.savefig("vnloads.jpg")
