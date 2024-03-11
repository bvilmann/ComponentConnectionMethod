# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 09:49:26 2024

@author: bvilm
"""
# ======================== PACKAGES ======================== 
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
# include parent directory
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# include ConnectionComponentMethod package
from CCM import ComponentConnectionMethod, StateSpaceSystem

# ======================== SIMULATION - CONTROL ======================== 
ccm_model = 'ccm_RLC'
file_path = f'{ccm_model}.xlsx'

# Parameters
R = .5  # Resistance in ohms
L = 0.01  # Inductance in henries
C = 0.001  # Capacitance in farads
Kp = 1
Ki = 30
params = {'R':R, 'L':L,'C':C,'Ki':Ki,'Kp':Kp}

# Load each subsystem and apply CCM
CCSM = [StateSpaceSystem(file_path, i, vals=params, comp_id=1) for i in (0,1)]
ccm = ComponentConnectionMethod(CCSM)

# ======================== SIMULATION - CONTROL ======================== 
# Time span for the simulation
t_span = (0, 0.5)  # 5 seconds
t_eval = np.linspace(*t_span, 5000)  # Time points where the solution is evaluated

# Input function (step input to 1V)
def input_function(t):
    return 1 if t > 0 else 0

# ------------------------ WITH PI CONTROL ------------------------
# Initial conditions
x0 = [0, 0, 0]  # Initial capacitor voltage and inductor current

# System of differential equations
def system_equations_ccm(t, x):
    u = input_function(t)
    return ccm.F @ x + ccm.G.flatten() * u

# Solve the system
solution_ccm = solve_ivp(system_equations_ccm, t_span, x0, t_eval=t_eval, method='RK45')

# ------------------------ NO CONTROL ------------------------
# System matrices
A = np.array([[0, 1/C],
              [-1/L, -R/L]])
B = np.array([[0],
              [1/L]])

# Initial conditions
x0 = [0, 0]  # Initial capacitor voltage and inductor current

# System of differential equations
def system_equations(t, x):
    u = input_function(t)
    return A @ x + B.flatten() * u

# Solve the system
solution = solve_ivp(system_equations, t_span, x0, t_eval=t_eval, method='RK45')

# ======================== PLOTTING ======================== 
# Plotting the results
fig, ax = plt.subplots(2,1,dpi=150,sharex=True)
ax[0].plot(solution.t, solution_ccm.y[1], label='With control',color='C0')
ax[1].plot(solution.t, solution_ccm.y[2], label='With control',color='C1')
ax[0].plot(solution.t, solution.y[0], label='No control',color='C0',ls=':')
ax[1].plot(solution.t, solution.y[1], label='No control',color='C1',ls=':')
ax[0].axhline(1,color='red',ls='--',label='Reference')
ax[1].set(xlabel='Time [s]',ylabel = 'Current [A]',xlim=t_span)
ax[0].set(ylabel = 'Voltage [V]')
for i in range(2):
    ax[i].grid(ls=':')
    ax[i].legend()
    ax[i].axhline(0, lw=0.75,color='k',ls='-')
fig.align_ylabels()
# plt.show()

plt.savefig('RLC_simulation.pdf',bbox_inches=0)
