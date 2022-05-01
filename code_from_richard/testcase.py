#!/usr/bin/env python3
"""
Numerical integration of a simpler set of equations with a constraint
to see if we get the same slow relaxation. We don't ––– but that's cos I
don't think the equations scale on the lattice... so the difference in rate
between the two processes is relevant here!
"""

import numpy as np, matplotlib.pyplot as plt, sys
from scipy.integrate import solve_ivp
from matplotlib.animation import FuncAnimation

M = 100 # Number of sites
until = 20
every = 0.1

system = 2

site1 = np.zeros(M-1)
site1[0] = 1

def odesys(t, y):
    # Create finite derivatives valid for n=1..M-1
    dyplus = (np.roll(y,-1) - y)[1:] # y[n+1]-y[n]
    dyminus = (np.roll(y,1) - y)[1:] # y[n-1]-y[n]

    if system == 0:
        return np.hstack(([0], dyplus+dyminus - y[1:] + site1)) # (0)
    elif system == 1:
        return np.hstack(([0], dyplus+dyminus + dyplus + np.sum(y) * site1)) # (1)
    return np.hstack(([0], dyplus+dyminus + dyplus / np.sum(y) + site1)) # (2)

# Initial condition should be [0] + [stuff]
# Needs to respect sum(y*np.arange(M)) = 1
y = np.zeros(M)
y[10] = 1/10 # Maybe delta function initial condition is

sol = solve_ivp(odesys, (0,until), y, t_eval=np.arange(0,until+every,every))

fig, ax = plt.subplots(1,1)
plot = ax.plot([], [], '-')[0]

txt = ax.text(0.5,0.95,'',ha='center',transform=ax.transAxes)
ax.set_xlim((0,M//3))
ax.set_ylim((0,np.max(sol.y)))
plt.xlabel(r'$k$')
plt.ylabel(r'$P(k)$')

ks = np.arange(M)

if system == 1:
    q = 0.5 # (1)
else:
    q = (3-np.sqrt(5))/2 # (0,2)

plt.plot(ks[1:], (1-q)**2 * q**(ks[1:]-1), '--')

def update(n):
    plot.set_data(ks[1:], sol.y[1:,n])
    txt.set_text(f"t={sol.t[n]:.2f}")
    return plot, txt

ani = FuncAnimation(fig, update, frames=len(sol.t), interval=50, repeat=True, blit=True, save_count=sys.maxsize)

plt.show()
