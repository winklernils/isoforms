#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 10:27:31 2023

This code runs the simulation of the full system model with two species
and nonlinear diffusion terms. The simulation is based on the FVM and the FiPy
package.

@author: Nils O. Winkler
"""

# Import necessary libraries
from fipy import Variable, CellVariable, Grid1D, TransientTerm, DiffusionTerm, \
    ExponentialConvectionTerm, ImplicitSourceTerm
from fipy.tools import numerix as nx
import numpy as np
from matplotlib import pyplot as plt
import os
import datetime
from methods import SAVEDATA, INTERPOLATE
import time

# Initialize simulation parameters
start_time = time.time()
BOOL_SAVE_RESULTS = True

# Saving parameters
save_every_step_lengths = 10
save_every_step_profiles = 1000
save_to_disk = 500
plot_every_step = 500

# Simulation constants
LENGTH_SIM = 1.0
MESH_nx = 50.0
MESH_dx = LENGTH_SIM / MESH_nx
TOTALTIME_SIM = 1.0
SWEEP_ACCURACY = 1e-8

# Generate simulation grid
mesh = Grid1D(nx=MESH_nx, dx=MESH_dx)
X_cell = mesh.cellCenters[0]
X = mesh.faceCenters[0]

# CFL-condition determines timestep and therefore accuracy
PERCENT_CFL_CONDITION = 50.0

# Define physical constants
PA = 0.05
PB = 0.18
K = 100.0
L_scr = 1.25
kr = 10.0
cr = 1.0
c_inf = 4.0

# Saving data
nonlinear_diffusion_type = "full_system"
COMMENT_pure = "paper_kymograph"
COMMENT_pure += "_every{}_{}th".format(save_every_step_lengths, 
                                       save_every_step_profiles) 

START_TIME = datetime.datetime.now()
TIMESTAMP = str(datetime.datetime.now()).split('.')[0]

HEADER = """# SIMULATION FULL MODEL ({})
# COEFFICIENTS: PA={}, PB={}, K={}, L_scr={}, C_inf={},
# NONLINEAR DIFF: TYPE: {},
# SIMULATION: N_mesh={}, Length_sim={}, Time_sim={}, SWEEP_ACC={}
""".format(TIMESTAMP, PA, PB, K, L_scr, c_inf, nonlinear_diffusion_type, 
           MESH_nx, LENGTH_SIM, TOTALTIME_SIM, SWEEP_ACCURACY)
CONSTANTS = [PA, PB, K, L_scr, c_inf, MESH_nx]

FILEPATH = "./data/{}/{}_k{}_pa{}_pb{}_L_scr{}_cinf{}".format(
    nonlinear_diffusion_type, COMMENT_pure, K, PA, PB, L_scr, c_inf)

# Create necessary filepaths
# Extract the directory path from FILEPATH
dir_path = os.path.dirname(FILEPATH)

# Create necessary directories if they do not exist
os.makedirs(dir_path, exist_ok=True)
print(f"Directories ensured for: {dir_path}")

if not os.path.exists(FILEPATH):
    os.mkdir(FILEPATH)
    print("Directory: " + FILEPATH +  " Created ")
else:    
    print("Directory: " + FILEPATH +  " already exists")


# Initial conditions
cA0 = 1.0
cB0 = cA0 / cr
L_0 = (1.0 + np.sqrt(1.0 - 4.0 * PA - 4.0 * PB / cr)) / 2.0 + 0.0001
G_0 = 0.0

# Scalar variables
L = Variable(name='L', value=L_0)  # cell length
L_plus = Variable(name='L_plus', value=G_0 + L_0 / 2)  # position of right cell edge
L_minus = Variable(name='L_minus', value=G_0 - L_0 / 2)  # position of left cell edge
G = Variable(name='G', value=G_0)  # position of cell's midpoint

# Old values for Euler stepping
L_old = L_0
L_plus_old = L_plus
L_minus_old = L_minus
G_old = G_0

# Field variables
sigma = CellVariable(name="stress", mesh=mesh, value=-L_0 * (L_0 - 1.0), 
                     hasOld=1)
cA = CellVariable(name="concentration", mesh=mesh, 
                  value=cA0 * (0.05 * (nx.cos(X_cell * nx.pi)) + 1.0), 
                  hasOld=1)  # cos shaped peak left
cB = CellVariable(name="concentration", mesh=mesh, value=cB0, hasOld=1)
L_plus_dot = K / L**2 * sigma.faceGrad[0, -1]
L_minus_dot = K / L**2 * sigma.faceGrad[0, 0]
L_dot = L_plus_dot - L_minus_dot
G_dot = (L_plus_dot + L_minus_dot) / 2.0

# Helper variables
v_tilde = -G_dot - L_dot * (X - 0.5) + 0.0 * sigma.faceGrad
u_dot = -K / (L**3) * sigma.faceGrad + L_dot / L * (X - 0.5) + G_dot / L

# Boundary conditions
sigma.constrain(-L * (L - 1.0), mesh.facesLeft)
sigma.constrain(-L * (L - 1.0), mesh.facesRight)
cA.faceGrad.constrain(0.0, where=mesh.exteriorFaces)
cB.faceGrad.constrain(0.0, where=mesh.exteriorFaces)

# Define the nonlinear diffusion term
Qtilde = (L * c_inf - (cA + cB))**2
DIFFUSION_TERM_CA_EQ = DiffusionTerm(coeff=1 / L * c_inf / Qtilde * 
                                     (L * c_inf - cB), var=cA)
DIFFUSION_TERM_CB_EQ = DiffusionTerm(coeff=kr / L * c_inf / Qtilde * 
                                     (L * c_inf - cA), var=cB)

# PDE system for stress and both concentrations
Qtilde = (L * c_inf - (cA.faceValue + cB.faceValue))**2

EQ_cA = TransientTerm(coeff=1., var=cA) == DIFFUSION_TERM_CA_EQ \
    - ExponentialConvectionTerm(coeff=1 / L * (K / L**2 * sigma.faceGrad + 
                                               v_tilde), var=cA) \
    + ExponentialConvectionTerm(coeff=1 / L * c_inf / Qtilde * cB.faceGrad, 
                                var=cA)

EQ_cB = TransientTerm(coeff=1., var=cB) == DIFFUSION_TERM_CB_EQ \
    - ExponentialConvectionTerm(coeff=1 / L * (K / L**2 * sigma.faceGrad + 
                                               v_tilde), var=cB) \
    + ExponentialConvectionTerm(coeff=kr / L * c_inf / Qtilde * cA.faceGrad, 
                                var=cB)

EQ_sigma = DiffusionTerm(coeff=L_scr / L**2, var=sigma) == \
    ImplicitSourceTerm(coeff=1., var=sigma) - (PA * cA + PB * cB)

eq = EQ_cA & EQ_cB & EQ_sigma

# Define lists for saving data
Ts = []
Ls = []
Gs = []
Gdots = []
L_dots = []

sigmas = []
concsA = []
concsB = []

escape_counter = 1  # initialize counter

# Save initial values
Ls.append(float(L.value))
L_dots.append(L_dot.value)
Gs.append(float(G.value))
Gdots.append(float(G_dot.value))
sigmas.append(np.array(sigma.value))
concsA.append(np.array((cA / L).value))
concsB.append(np.array((cB / L).value))

# Courant-Friedrichs-Lewy (CFL) condition to determine timestep dt
MAX_VEL_HAT = np.max([np.max(np.abs(u_dot.value * L.value)), 0])
diffA = max(1 / L.value * c_inf / Qtilde * (L.value * c_inf - cB.faceValue))
diffB = max(kr / L.value * c_inf / Qtilde * (L.value * c_inf - cA.faceValue))
MAX_DIFF_C_NONLIN = max(diffA, diffB)
MAX_DIFF_C_NONLIN = 1.0
dt = PERCENT_CFL_CONDITION * MESH_dx**2 * L.value**2 / (MAX_VEL_HAT * 
        L.value * MESH_dx + 2.0 * np.max([MAX_DIFF_C_NONLIN, np.abs(L_scr)]))

# Initialize total time
total_t = 0.0
Ts.append(total_t)

# Simulation loop
while total_t < TOTALTIME_SIM:
    """
    The algorithm goes as follows with the values at timepoint (n) and (n+1):
    0) We start with L(n), sigma(n), cA(n), cB(n) as ICs and L_dot(n), G_dot(n) as function of them
    1) Calculate L(n+1) = L(n) + L_dot(n)*dt and G(n+1) as Euler step
    2) Solve PDEs for cA(n+1), cB(n+1), sigma(n+1), L_dot(n+1), G_dot(n+1)
        with L(n+1), i.e. BC fulfilled.
        Consider it is a implicit equation with sigma(n+1) on both sides --> sweeping
    3) repeat
    """
    sigma.updateOld()
    cA.updateOld()
    cB.updateOld()
    
    # Calculate new length L(n+1)
    if total_t == 0.0:
        L.setValue(L_old)
        L_plus.setValue(L_plus_old)
        L_minus.setValue(L_minus_old)
        G.setValue(G_old)
    else:
        L.setValue(L_old + L_dot.value * dt)
        G.setValue(G_old + G_dot.value * dt)
        L_plus.setValue(L_plus_old + L_plus_dot.value * dt)  # separate data for edges to plot phase relation
        L_minus.setValue(L_minus_old + L_minus_dot.value * dt)

    L_old = L.value
    G_old = G.value
    L_plus_old = L_plus.value
    L_minus_old = L_minus.value
    
    # Sweep the coupled equation until residual is small enough for sigma(n+1), cA(n+1), cB(n+1)
    RESIDUAL = 100.0
    sweep_counter__ = 0
    while RESIDUAL >= SWEEP_ACCURACY and sweep_counter__ < 100:
        RESIDUAL = eq.sweep(dt=dt)
        sweep_counter__ += 1
        if sweep_counter__ == 99:
            print("WARNING! RESIDUALS NOT CONVERGING after " + 
                  str(sweep_counter__) + " sweeps at timestep: " + 
                  str(escape_counter) + ". Res.: " + str(RESIDUAL))
    
    total_t += dt
    
    # Save values
    if escape_counter % save_every_step_lengths == 1 or save_every_step_lengths == 1:
        Ts.append(total_t)
        Ls.append(float(L.value))
        L_dots.append(L_dot.value)
        Gs.append(float(G.value))
        Gdots.append(float(G_dot.value))
        
    if escape_counter % save_every_step_profiles == 1 or save_every_step_profiles == 1:
        sigmas.append(np.array(sigma.value))
        concsA.append(np.array((cA / L).value))
        concsB.append(np.array((cB / L).value))
        
    # Plotting examples of the fields
    if escape_counter % plot_every_step == 1 or plot_every_step == 1:
        fig = plt.figure(figsize=(6, 6))
        fig.add_subplot(221)
        plt.plot(X_cell.value, cA.value, label=r"$LcA$")
        plt.plot(X_cell.value, cB.value, label=r"$LcB$")
        plt.legend()
        fig.add_subplot(222)
        plt.plot(Ts, np.asarray(Ls), label="L")
        plt.legend()
        fig.add_subplot(223)
        plt.plot(Ts, Gdots, label='Gdot')
        plt.legend()
        fig.add_subplot(224)
        plt.plot(X_cell.value, sigmas[-1], label=r'$\sigma$')
        plt.legend()
        
        plt.show()
        print("sweep counter = " + str(sweep_counter__))
        print("total time = " + str(total_t))
        print("steps so far = " + str(escape_counter))
        print("-------------------------------")
        
    # Save to disk every 5000 timesteps
    if escape_counter % save_to_disk == 0:
        datatosave = [np.array(X_cell.value), np.array(sigmas), 
                      np.array(concsA), np.array(concsB), np.array(Ts), 
                      np.array(Ls), np.array(L_dots), np.array(Gs), 
                      np.array(Gdots), np.array(X.value)]
        namestosave = ["X", "sigma", "cA", "cB", "T", "L", "Ldots", "G", 
                       "Gdots", "X_face"]
        if BOOL_SAVE_RESULTS: SAVEDATA(datatosave, namestosave, FILEPATH, 
                                       HEADER)
        
    # Update CFL condition
    MAX_VEL_HAT = np.max([np.max(np.abs(u_dot.value * L.value)), 0])
    diffA = max(1 / L.value * c_inf / Qtilde * (L.value * c_inf - 
                                                cB.faceValue))
    diffB = max(kr / L.value * c_inf / Qtilde * (L.value * c_inf - 
                                                 cA.faceValue))
    MAX_DIFF_C_NONLIN = max(diffA, diffB)
    MAX_DIFF_C_NONLIN = 1.0
    dt_cfl = PERCENT_CFL_CONDITION * MESH_dx**2 * L.value**2 / (MAX_VEL_HAT * 
                L.value * MESH_dx + 2.0 * np.max([MAX_DIFF_C_NONLIN, np.abs(L_scr)]))
    dt = dt_cfl
    
    escape_counter += 1

# Save the final data to disk
datatosave = [np.array(X_cell.value), np.array(sigmas), np.array(concsA), 
              np.array(concsB), np.array(Ts), np.array(Ls), np.array(L_dots), 
              np.array(Gs), np.array(Gdots), np.array(X.value)]
namestosave = ["X", "sigma", "cA", "cB", "T", "L", "Ldots", "G", "Gdots", 
               "X_face"]
if BOOL_SAVE_RESULTS: SAVEDATA(datatosave, namestosave, FILEPATH, HEADER)

end_time = time.time()
total_time = end_time - start_time
print(total_time)

plt.plot(Ts, Gdots)
print(len(Ls))
print(Ls[-200:])