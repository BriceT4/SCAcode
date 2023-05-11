#!/usr/bin/python
# SCAcode_Setup.py
# Base of SCAcode project
# (C) Brice Turner, 2023

import argparse
import importlib
import importlib.util
import os

# BEGIN: COMMAND LINE INTERFACE ##############################################
def SetupCommandLine():
    print(f"""
##############################################################################
####       #######      ####       ###########################################
##    ###   ####    #######   ###   #############################  ###########
#    ##########   ########   #####   ############################  ###########
#    #########   #########  #######  ############################  ###########
##    ########   #########  #######  ############################  ###########
####     #####   #########           ###       ###      ####       ###      ##
######    ####   #########    ###    ##   #######   ##   ##   ###  ##   ##   #
######    #####   ########   #####   ##  ########  ####  ##  ####  ##   ##   #
#   #    #######    ######   #####   ##   #######   ##   ##   ###  ###     ###
##     ###########      ##   #####   ###       ###      ####       ####      #
##############################################################################
""")

    desc = f'Command line interface for {os.path.basename(__file__)}.'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--INPUT', nargs='?', dest='infile',
                        metavar='input_file', required=True,
                        help='Input file name. REQUIRED. Must be .py.\n'
                            'Usage: -i/--INPUT <filename>')
    args = parser.parse_args()

    # Get the file path and the file name without the extension
    file_path = os.path.abspath(args.infile)
    infile_name = os.path.splitext(os.path.basename(args.infile))[0]

    # Replace periods/dots with underscores in the module name
    infile_name = infile_name.replace('.', '_')

    # Make sure the input file exists
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"File '{args.infile}' not found")

    # Import the input file as a module
    spec = importlib.util.spec_from_file_location(infile_name, file_path)
    input_file = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(input_file)

    return input_file, infile_name
# END:   COMMAND LINE INTERFACE ##############################################



# BEGIN: SETUP ###############################################################
def solveIt(input_file, infile_name):
    from src.SCAcode import solver
    fp_data, dir_output = solver(input_file, infile_name)
    return fp_data, dir_output
# END:   SETUP ###############################################################



# BEGIN: PLOT ################################################################
def plotIt(fp_data, infile_name, dir_output):
    from src.genPlots import plotter
    plotter(fp_data, infile_name, dir_output)
    return
# END:   PLOT ################################################################



# BEGIN: SUCCESS MESSAGE #####################################################
def PrintSuccess():
    print(f"""
##############################################################################
{os.path.basename(__file__)} sucessfully ran.
##############################################################################
    """)
# END:   SUCCESS MESSAGE #####################################################



# BEGIN: EXECUTE #############################################################
input_file, infile_name = SetupCommandLine()
fp_data, dir_output = solveIt(input_file, infile_name)
plotIt(fp_data, infile_name, dir_output)
PrintSuccess()
# END:   EXECUTE #############################################################
