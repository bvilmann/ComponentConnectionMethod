# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 09:49:26 2024
@author: bevil @ Technical University of Denmark (DTU)
"""
# ======================== PACKAGES AND FUNCTIONS ========================
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
# include parent directory
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# include ConnectionComponentMethod package
from CCM import ComponentConnectionMethod, StateSpaceSystem
# System of differential equations
def system_equations_ccm(t, x):
    u = input_function(t)
    return ccm.F @ x + ccm.G.flatten() * u

# Input function (step input to 1V)
def input_function(t):
    return 1 if t > 0 else 0
    
    

#%% ======================== TIME DOMAIN COMPARISON ======================== 
models = ['ccm_RLC','ccm_RLC_PIR']
# Parameters
R = 1  # Resistance in ohms
L = 0.01  # Inductance in henries
C = 0.001  # Capacitance in farads
Kp = 1
Ki = 30
Kr = 1
w_c = 1400
params = {'R':R, 'L':L,'C':C,'Ki':Ki,'Kp':Kp,'Kr':Kr,'w_c':w_c}
x0 = [0, 0, 0]  # Initial capacitor voltage and inductor current
t_span = (0, 0.5)  # 5 seconds
t_eval = np.linspace(*t_span, 5000)  # Time points where the solution is evaluated

fig, ax = plt.subplots(2,1,dpi=150,sharex=True)
for i, model in enumerate(models):
    x0 = [0,0,0] if i == 0 else [0,0,0,0]
    # Load each subsystem and apply CCM
    file_path = f'{os.path.dirname(os.path.abspath(__file__))}\\{model}.xlsx'
    ccm = ComponentConnectionMethod(file_path,params=params)

    print(ccm)

    # Solve the system
    solution_ccm = solve_ivp(system_equations_ccm, t_span, x0, t_eval=t_eval, method='RK45')

    ax[0].plot(t_eval, solution_ccm.y[1], label=model,color=['C0','C4'][i],ls=['-','-'][i])
    ax[1].plot(t_eval, solution_ccm.y[2], label=model,color=['C0','C4'][i],ls=['-','-'][i])

    del ccm

ax[1].set(xlabel='Time [s]',ylabel = 'Current [A]',xlim=t_span)
ax[0].set(ylabel = 'Voltage [V]')
ax[0].axhline(1,color='k',ls='--',lw=0.75,label='Reference')
for i in range(2):
    ax[i].grid(ls=':')
    ax[i].legend()
    ax[i].axhline(0, lw=0.75,color='k',ls=':')
fig.align_ylabels(ax)

fig.tight_layout()

plt.show()



#%% ======================== EIGENVALUE SENSITIVY ANALYSIS ======================== 
fig, ax1 = plt.subplots(1,1,dpi=150)

file_path = f'{os.path.dirname(os.path.abspath(__file__))}\\CCM_RLC.xlsx'

for i, Ki in enumerate([0,10,30,50]):
    
    # Dynamically define the parameters of interest
    params = {'R':R, 'L':L,'C':C,'Ki':Ki,'Kp':Kp}
    
    # Load each subsystem and apply CCM
    ccm = ComponentConnectionMethod(file_path,params=params)
    print(ccm)
    
    ccm.show(fontsize = 24,save=True,boldfacecolor='magenta',text_alpha = 0.75)
        
    eigs = np.linalg.eigvals(ccm.F)
    ax1.scatter(eigs.real,eigs.imag,color='C0',label=Ki,alpha=1/(i+1))

ax1.legend(title='$K_i$')
ax1.grid(ls=':')
plt.show()


    