#!/usr/bin/python
# <input>.py
# An input file for SCAcode project
# (C) 2023, Brice Turner and Parker McLeod

L = 3.0 # (m)
L_e = 3.45 # (m)
num_CV = 400 # number of control volumes along length. 

D = 10/1000 # (m), converted from (mm)
D_ci = 9.1/1000 # (m), converted from (mm) 
D_fo = 8.9/1000 # (m), converted from (mm) 
Pitch = 12.9/1000 # (m), converted from (mm)

T_m_in = 260 # (Celcius)
P_nom = 6900000 # (Pa)
qp_max = 30*1000 # (W/m) converted from (kW/m) 
m_dot = 0.17 # (kg/s)
k_cl = 16 # cladding conductivity (W/m-K)

g = 9.81 # gravitational acceleration on Earth Z(m/s^2)
SB_constant = 5.67E-8 # Stefan-Boltzmann constant (W m^-2 K^-4)

Reactor_Type = 'BWR'
SCB_flag = True
CHFR_crit_limit = 1.9

convergence_crit_gap = 0.001
convergence_crit_fuel = 0.0001 

