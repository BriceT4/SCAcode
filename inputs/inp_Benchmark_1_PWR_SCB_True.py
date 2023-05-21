#!/usr/bin/python
# <input>.py
# An input file for SCAcode project
# (C) 2023, Brice Turner and Parker McLeod

L = 3.7 # (m)
L_e = 4 # (m)
num_CV = 400 # number of control volumes along length. 

D = 9.4/1000 # (m), converted from (mm)
D_ci = 8.5/1000 # (m), converted from (mm) 
D_fo = 8.4/1000 # (m), converted from (mm) 
Pitch = 12.5/1000 # (m), converted from (mm)

T_m_in = 280 # (Celcius)
P_nom = 15600000 # (Pa)
qp_max = 37*1000 # (W/m) converted from (kW/m) 
m_dot = 0.30 # (kg/s)
k_cl = 16 # cladding conductivity (W/m-K)

g = 9.81 # gravitational acceleration on Earth Z(m/s^2)
SB_constant = 5.67E-8 # Stefan-Boltzmann constant (W m^-2 K^-4)

Reactor_Type = 'PWR'
SCB_flag = True
CHFR_crit_limit = 1.1

convergence_crit_gap = 0.001
convergence_crit_fuel = 0.0001 

