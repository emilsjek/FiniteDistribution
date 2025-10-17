# -*- coding: utf-8 -*-
"""
Author: Emils Jekabsons
"""

import numpy as np
import matplotlib.pyplot as plt

#Given parameters
L = 0.2           
T0 = 475
TL = 290          
Ts = 290          
hc = 40           
P = 0.015         
k = 240           
Ac = 1.55e-5      

#Derived coefficient
alpha = (hc * P) / (k * Ac)

#Discretization
n = 20 # number of segments
dx = L / n
x = np.linspace(0, L, n+1)

#Matrix assembly
A = np.zeros((n-1, n-1))
b = np.zeros(n-1)

for i in range(n-1):
    if i != 0:
        A[i, i-1] = 1
    A[i, i] = -2 - alpha * dx**2
    if i != n-2:
        A[i, i+1] = 1
    b[i] = -alpha * dx**2 * Ts  # RHS term

# Apply boundary conditions
b[0] -= T0
b[-1] -= TL

#Solve for internal nodes
T_internal = np.linalg.solve(A, b)

# Combine with boundary nodes
T = np.zeros(n+1)
T[0] = T0
T[-1] = TL
T[1:-1] = T_internal

#Plot results
plt.figure(figsize=(8,5))
plt.plot(x, T, 'r-', marker='o', label='Temperature Distribution')
plt.xlabel('Fin Length (m)')
plt.ylabel('Temperature (K)')
plt.title('Temperature Distribution Along the Pin Fin')
plt.grid(True)
plt.legend()
plt.show()

#Output data
for i in range(len(x)):
    print(f"x = {x[i]:.3f} m,  T = {T[i]:.2f} K")