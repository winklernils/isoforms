#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 10:27:31 2023

This script simulates the model with two species and a rigid limit
including nonlinear diffusion. The model is solved with the finite volume method fipy.

@author: nilswinkler
"""

# Import necessary libraries
from fipy import Variable, CellVariable, Grid1D, TransientTerm, DiffusionTerm, \
    UpwindConvectionTerm, ImplicitSourceTerm
from fipy.tools import numerix as nx
import numpy as np
from matplotlib import pyplot as plt
import os
import datetime
from methods import SAVEDATA, INTERPOLATE
import scipy

# Define physical constants
PA = 0.05
PB = 0.18
chir = PA / PB
L_scr = 1.25
cr = 1.0
kr = 4.0
lam = 44.0
c_inf = 4.0

BOOL_SAVE_RESULTS = True

# Saving parameters
save_every_step_lengths = 5
save_every_step_profiles = 5
save_to_disk = 1000
plot_every_step = 500


# Simulation constants
LENGTH_SIM = 1.0
MESH_nx = 50.0
MESH_dx = LENGTH_SIM / MESH_nx
TOTALTIME_SIM = 1.8
SWEEP_ACCURACY = 1e-6

# Generate simulation grid
mesh = Grid1D(nx=MESH_nx, dx=MESH_dx)
X_cell = mesh.cellCenters[0]
X = mesh.faceCenters[0]

# CFL-condition determines timestep and therefore accuracy
PERCENT_CFL_CONDITION = 50.0

# Saving data
nonlinear_diffusion_type = 'rigid_limit'
COMMENT_pure = "paper"
COMMENT_pure += "_every{}_{}th".format(save_every_step_lengths, 
                                       save_every_step_profiles) 

START_TIME = datetime.datetime.now()
TIMESTAMP = str(datetime.datetime.now()).split('.')[0]

HEADER = """# SIMULATION RECHO MODEL ({})
# COEFFICIENTS: lam={}, chir={}, kr={}, C_inf={},
# NONLINEAR DIFF: TYPE: {},
# SIMULATION: N_mesh={}, Length_sim={}, Time_sim={}, SWEEP_ACC={}
""".format(TIMESTAMP, lam, chir, kr, c_inf, nonlinear_diffusion_type, 
           MESH_nx, LENGTH_SIM, TOTALTIME_SIM, SWEEP_ACCURACY)
CONSTANTS = [chir, lam, kr, c_inf, MESH_nx]

FILEPATH = "./data/rigid_L/{}_K{}_kr{}_pa{}_pb{}_L_scr{}_cinf{}".format(
    COMMENT_pure, lam, kr, PA, PB, L_scr, c_inf)

# Create necessary filepaths
dir_path = os.path.dirname(FILEPATH)
os.makedirs(dir_path, exist_ok=True)
print(f"Directories ensured for: {dir_path}")

if not os.path.exists(FILEPATH):
    os.mkdir(FILEPATH)
    print("Directory: " + FILEPATH + " Created ")
else:
    print("Directory: " + FILEPATH + " already exists")

# Initial conditions
cA0 = 1.0
cB0 = cA0 / cr
s0 = lam * (1.0 + 1.0 / (chir * cr))
G_0 = 0.0

# Scalar variables
G = Variable(name='G', value=G_0)  # position of cell's midpoint
G_old = G_0  # old value for Euler stepping

# Field variables
s = CellVariable(name="stress", mesh=mesh, value=0.05, hasOld=1)  # small nonzero value needed to get dynamics going
cA = CellVariable(name="concentration", mesh=mesh, 
                  value=cA0 * (0.05 * (nx.cos(X_cell * nx.pi)) + 1.0), 
                  hasOld=1)  # cos shaped peak left
cB = CellVariable(name="concentration", mesh=mesh, value=cB0, hasOld=1)
Gdot = s.faceGrad[0, 0]

# Helper variables
u_dot = Gdot - s.faceGrad

# Boundary conditions
s.constrain(0.0, mesh.facesLeft)
s.constrain(0.0, mesh.facesRight)
cA.faceGrad.constrain(0.0, where=mesh.exteriorFaces)
cB.faceGrad.constrain(0.0, where=mesh.exteriorFaces)

# Define the nonlinear diffusion term
Q = (c_inf - (cA + cB))**2
DIFFUSION_TERM_CA_EQ = DiffusionTerm(coeff=c_inf / Q * (c_inf - cB), var=cA)
DIFFUSION_TERM_CB_EQ = DiffusionTerm(coeff=kr * c_inf / Q * (c_inf - cA), 
                                     var=cB)

# PDE system for stress and both concentrations
Q = (c_inf - (cA.faceValue + cB.faceValue))**2

EQ_cA = TransientTerm(coeff=1., var=cA) == DIFFUSION_TERM_CA_EQ \
    + UpwindConvectionTerm(coeff=(Gdot - s.faceGrad), var=cA) \
    + UpwindConvectionTerm(coeff=c_inf / Q * cB.faceGrad, var=cA)

EQ_cB = TransientTerm(coeff=1., var=cB) == DIFFUSION_TERM_CB_EQ \
    + UpwindConvectionTerm(coeff=Gdot - s.faceGrad, var=cB) \
    + UpwindConvectionTerm(coeff=kr * c_inf / Q * cA.faceGrad, var=cB)

EQ_s = DiffusionTerm(coeff=L_scr, var=s) == ImplicitSourceTerm(coeff=1., var=s) \
    + s0 - lam * (cA + cB / chir)

eq = EQ_cA & EQ_cB

# Define lists for saving data
Ts = []
Gs = []
Gdots = []
u_dots = []
s0s = []
ss = []
concsA = []
concsB = []

escape_counter = 1  # initialize counter

# Save initial values
Gs.append(float(G.value))
Gdots.append(float(0))  # initial zero velocity
s0s.append(s0)
concsA.append(np.array(cA.value))
concsB.append(np.array(cB.value))

# Courant-Friedrichs-Lewy (CFL) condition to determine timestep dt
MAX_VEL_HAT = np.max([np.max(np.abs(u_dot.value)), 0])
diffA = max(c_inf / Q * (c_inf - cB.faceValue))
diffB = max(kr * c_inf / Q * (c_inf - cA.faceValue))
MAX_DIFF_C_NONLIN = max(diffA, diffB)
dt = PERCENT_CFL_CONDITION * MESH_dx**2 / (MAX_VEL_HAT * MESH_dx 
        + 2.0 * np.max([MAX_DIFF_C_NONLIN, np.abs(L_scr)]))
total_t = 0.0
Ts.append(total_t)

def solve_stress_eq(s0, ca, cb, s):
    """Solve the stress equation temporarily with given s0, boundary conditions
    and concentrations. Needed because s0 changes with every time step
    without analytic solution. 
    Output: stress gradient at left and right boundary, stress field
    
    Difference between left and right subject to root finding later on 
    """
    s0 = s0[0] if isinstance(s0, np.ndarray) else s0
    s0_var = Variable(name='s0', value=s0)
    s_temp = CellVariable(name="stress", mesh=mesh, value=s.value, hasOld=1)
    s_temp.constrain(0.0, mesh.facesLeft)
    s_temp.constrain(0.0, mesh.facesRight)
    EQ_s = DiffusionTerm(coeff=L_scr, var=s_temp) == ImplicitSourceTerm(coeff=1., 
            var=s_temp) + s0_var - lam * (ca + cb / chir)
    EQ_s.solve(var=s_temp)
        
    spl = s_temp.faceGrad.value[0, 0]  # left side
    spr = s_temp.faceGrad.value[0, -1]  # right side
    
    return (spl - spr, s_temp)

no_convergence_counter = 0    
VALUE_CONVERGENCE = False
DIVERGING_RESIDUALS = False

# Simulation loop
while total_t < TOTALTIME_SIM:
    """
    The algorithm goes as follows with the values at timepoint (n) and (n+1):
    0) We start with G(n), sigma(n), cA(n), cB(n) as ICs and Gdot(n) as function of them
    1) determine s0(n) that satisfies BCs by solving the stress equation with 
        cA(n), cB(n), sigma(n)
    2) Calculate G(n+1) as Euler step
    3) Solve PDEs for cA(n+1), cB(n+1), sigma(n+1), G_dot(n+1)
    4) repeat
    """
    s.updateOld()
    cA.updateOld()
    cB.updateOld()
    
    # Adapt s0 for matching stress derivatives at boundary
    eq_temp = lambda s0: solve_stress_eq(s0, cA, cB, s)[0]
    res = scipy.optimize.root_scalar(eq_temp, x0=s0)  # find root of eq_temp   
    s0 = res.root  # overwrite s0 with better fit at found root
    
    s_new = solve_stress_eq(s0, cA, cB, s)[1]  # corrected stress field
    s.setValue(s_new)  # overwrite stress field with corrected value

    # Calculate new velocity G(n+1)
    if total_t == 0.0:
        G.setValue(G_old)
    else:
        G.setValue(G_old + Gdot.value * dt)
        
    G_old = G.value
    
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
            no_convergence_counter += 1
            
    if no_convergence_counter > 10:
        print("WARNING! RESIDUALS NOT CONVERGING after 10 times. BREAKING SIMULATION")
        DIVERGING_RESIDUALS = True
        break
    total_t += dt
    
    # Save values
    if escape_counter % save_every_step_lengths == 1 or save_every_step_lengths == 1:
        Ts.append(total_t)
        Gs.append(float(G.value))
        Gdots.append(float(Gdot.value))
        s0s.append(s0)
        
    if escape_counter % save_every_step_profiles == 1 or save_every_step_profiles == 1:
        concsA.append(np.array(cA.value))
        concsB.append(np.array(cB.value))
        ss.append(np.array(s.value))
    
    # Plotting examples of the fields
    if escape_counter % plot_every_step == 1 or plot_every_step == 1:
        fig = plt.figure(figsize=(6, 6))
        fig.add_subplot(221)
        plt.plot(X_cell.value, cA.value, label=r"$LcA$")
        plt.plot(X_cell.value, cB.value, label=r"$LcB$")
        plt.legend()
        fig.add_subplot(222)
        plt.plot(Ts, np.asarray(s0s), label=r"$s_0$")
        plt.legend()
        fig.add_subplot(223)
        plt.plot(Ts, Gdots, label='Gdot')
        plt.legend()
        fig.add_subplot(224)
        plt.plot(X_cell.value, s.value, label=r'$s$')
        plt.legend()
        
        plt.show()
        print("sweep counter = " + str(sweep_counter__))
        print("total time = " + str(total_t))
        print("steps so far = " + str(escape_counter))
        print("-------------------------------")
    
    
        
    # Update CFL condition
    MAX_VEL_HAT = np.max([np.max(np.abs(u_dot.value)), 0])
    diffA = max(c_inf / Q * (c_inf - cB.faceValue))
    diffB = max(kr * c_inf / Q * (c_inf - cA.faceValue))
    MAX_DIFF_C_NONLIN = max(diffA, diffB)
    MAX_DIFF_C_NONLIN = 1.0
    
    dt_cfl = PERCENT_CFL_CONDITION * MESH_dx**2 / (MAX_VEL_HAT * MESH_dx 
                + 2.0 * np.max([MAX_DIFF_C_NONLIN, np.abs(L_scr)]))
    if dt_cfl <= dt: 
        dt = dt_cfl
    
    escape_counter += 1

# Save the final data to disk
datatosave = [np.array(s0s), np.array(Ts), np.array(Gs), np.array(Gdots), 
              np.array(concsA), np.array(concsB), np.array(X.value)]
namestosave = ["s0", "T", "G", "Gdots", "cA", "cB", "X_face"]
if BOOL_SAVE_RESULTS: SAVEDATA(datatosave, namestosave, FILEPATH, HEADER)
