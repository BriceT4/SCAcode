#!/usr/bin/python
# genPlots.py
# Plotter of of SCAcode project
# (C) Brice Turner, 2023

import matplotlib.pyplot as plt
import os
import pandas as pd
import time

time_start = time.time()

def plotter(fp_data, infile_name, dir_output):

    data = pd.read_csv(fp_data)
    data = data.to_numpy()

    z = data[0:-1, 0]
    T_m = data[0:-1, 1]
    T_co = data[0:-1, 2]
    T_ci = data[0:-1, 3]
    T_fo = data[0:-1, 4]
    T_max = data[0:-1, 5]
    x = data[0:-1, 6]
    x_e = data[0:-1, 7]
    CHFR = data[0:-1, 10]
    Boiling_out_flag = data[0:-1, 9]
    SCB_out_flag = data[0:-1, 10]
    DeltaP_thisCell = data[0:-1, 12]
    Re = data[0:-1, 13]



    fig,ax = plt.subplots(3, figsize=(10,8))
    ax[0].scatter(z, T_max, marker = '.', facecolors='k', edgecolors='k', label = '$T_{max}$')
    ax[0].scatter(z, T_fo, marker = 's', facecolors='m', edgecolors='m', label = '$T_{fo}$')
    ax[0].set_ylabel('Temperature ($^{\circ}$C)')
    ax[0].legend(loc = 'upper left')
    # ax[0].legend(loc = [0.9, 0.6])

    ax[1].scatter(z, T_ci, marker = '.', facecolors='k', label = '$T_{ci}$')
    ax[1].scatter(z, T_co, marker = 's', facecolors='m', label = '$T_{co}$')
    ax[1].scatter(z, T_m, marker = '^', facecolors='b', label = '$T_m$')
    ax[1].set_ylabel('Temperature ($^{\circ}$C)')
    # ax[1].legend(loc = 'upper left')
    ax[1].tick_params(labelright=True)
    ax1 = ax[1].twinx()
    ax1.scatter(z, CHFR, marker = '+', facecolors='r', label = 'CHFR')
    ax1.set_ylabel('CHFR (arb. unit)', color='r')
    ax1.tick_params(axis='y', colors='r')
    h1, l1 = ax[1].get_legend_handles_labels()
    h2, l2 = ax1.get_legend_handles_labels()
    ax1.legend(h1+h2, l1+l2, loc = 'upper left')
    ax1.tick_params(color = 'r')

    ax[2].scatter(z, x, marker = '.', facecolors='k', label = '$x$')
    ax[2].scatter(z, x_e, marker = 's', facecolors='m', label = '$x_e$')
    ax[2].set_ylabel('Quality (arb. unit)')
    ax[2].set_xlabel('z (m)')
    ax2 = ax[2].twinx()
    ax2.scatter(z, Re, marker = '^', facecolors='r', label = 'Re')
    ax2.set_ylabel('Reynolds nubmer (arb. unit)', color='r')
    ax2.tick_params(axis='y', colors='r')
    h1, l1 = ax[2].get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax2.legend(h1+h2, l1+l2, loc = 'lower right')
    ax2.tick_params(color = 'r')


    plt.tight_layout()
    fig = plt.gcf()

    timestr = time.strftime('%Y%m%d_%H%M%S')
    filename = f'SCAcode_plots_{infile_name}_{timestr}.png'

    if not os.path.exists(dir_output):
        os.makedirs(dir_output)
    fp_data = os.path.join(dir_output, filename)
    fig.savefig(fp_data, dpi=600)

    time_end = time.time()
    time_elapsed = time_end - time_start

    print(f"""
{os.path.basename(__file__)} sucessfully ran.
Elapsed run time: {time_elapsed:.2f} s.
Results were exported to {fp_data}.
    """)

    # plt.show()
    return

