#!/usr/bin/python
# SCAcode.py
# Equation solver of SCAcode project
# (C) 2023, Brice Turner and Parker McLeod

import csv
import numpy as np
import os
import time

from scipy.interpolate import interp1d

# BEGIN: SOLVER FUNCTION #####################################################
def solver(inp, infile_name):
    time_start = time.time()

    def read_data(file_path):
        with open(file_path, 'r') as file:
            csv_reader = csv.reader(file, delimiter='\t')
            header = ['T (C)', 'P_sat (Pa)', 'vol_l (m^3/kg)', 'vol_g (m^3/kg)',
                    'h_l (J/kg)', 'h_g (J/kg)', 'mu_l (kg/m-s)', 'k_l (W/m-K)',
                    'Pr_l (arb. unit)', 'mu_g (kg/m-s)']
            proptable = [dict(zip(header, map(float, row))) for row in csv_reader]  

        return proptable
    proptable = read_data('.\\data\\proptable.txt')

    # BEGIN: COMPUTE PROPERTY INTERPOLATIONS #################################
    def interpPropTable(YourProp, YourValue, proptable):
        lower_row = None
        upper_row = None
        for row in proptable:
            if row[YourProp] == YourValue:
                return row
            elif row[YourProp] < YourValue:
                lower_row = row
            else:
                upper_row = row
                break
        if lower_row is None or upper_row is None:

            return None
        interp_values = {prop: interp1d([lower_row[YourProp], upper_row[YourProp]],
                                        [lower_row[prop], upper_row[prop]])(YourValue) for prop in proptable[0].keys()}
        interp_values['rho_l (kg/m^3)'] = 1 / interp_values['vol_l (m^3/kg)']
        interp_values['rho_g (kg/m^3)'] = 1 / interp_values['vol_g (m^3/kg)']

        return interp_values
        # END:   COMPUTE PROPERTY INTERPOLATIONS #################################

    # BEGIN: VARIOUS CALCULATIONS AND VARIABLE INITIALIZATIONS ###############
    # define saturation properties needed later,
    # based on a constant pressure from input file
    data = []
    props_sat = interpPropTable('P_sat (Pa)', inp.P_nom, proptable) 
    T_sat = props_sat['T (C)']
    h_l_sat = props_sat['h_l (J/kg)']
    h_g_sat = props_sat['h_g (J/kg)']

    # create discretized length vector
    Delta_z = inp.L/inp.num_CV
    z_vector = np.linspace(-inp.L/2 + Delta_z, inp.L/2, inp.num_CV)

    #calculate relevant input-file-based constants of the problem.
    R_co = inp.D/2 # radius, cladding outer (m)
    R_ci = inp.D_ci/2 # radius, cladding inner (m)
    R_fo = inp.D_fo/2 # radius, fuel inner (m)
    R_g = (R_fo + R_ci)/2 # radius, middle of fuel-cladding gap (m)
    delta_eff = R_ci - R_fo # effective cladding-fuel gap (m)
    Area = inp.Pitch**2 - np.pi/4*inp.D**2 # cross-sectional area of interest (m^2)
    D_e = inp.D*(4/np.pi*(inp.Pitch/inp.D)**2 - 1) # effective/hydraulic diameter (m)
    G = inp.m_dot/Area # mass flux (kg/m^2-s)

    #initalize some boiling variables
    x = 0 # quality
    z_D = 0 # first z locations where T_co > T_sat 
    x_e_zD = 0 # equilibrium quality at z_D
    CHFR_crit_flag = False # "CHFR has been reached" flag
    first_zD_flag = True # "z_D has not been reached flag"

    # initialize h_i+1 to <not a number> so once
    # the calculations using it starts, it will seemlessly be 
    # converted to a number and not cause issues in the first CV
    h_iplus1 = False 
    # END:   VARIOUS CALCULATIONS AND VARIABLE INITIALIZATIONS ###############



    # BEGIN: DEFINE q' CALCULATION FUNCTIONS #################################
    # define function to calculate q'(z)_{i+1/2} as a function of z.
    def qp_iPlusHalf(YourZ):
        qp_iplushalf = inp.qp_max*np.cos(np.pi*(YourZ - 0.5*Delta_z)/inp.L_e)

        return(qp_iplushalf)

    # define function to calculate q\prime(z)_{i+1} as a function of z
    def qp_iPlus1(YourZ):
        qp_iplus1 = inp.qp_max*np.cos(np.pi*(YourZ)/inp.L_e)

        return(qp_iplus1)
    # END:   DEFINE q' CALCULATION FUNCTIONS #################################



    # BEGIN: OUTER LOOP OF Z VECTOR ##########################################S
    for z in z_vector:
        print(f'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% z =  {z}\n')

        # BEGIN: COMPUTE ENERGY BALANCE IN COOLANT Tm(z) #####################
        # if you're at the very bottom of the channel,
        # use values from input file.
        if z == -inp.L/2 + Delta_z:
            # your temp is the input temp
            T_m = inp.T_m_in
        # we can now move out of the "if:" since it 
        # doesn't matter where we are in the channel, 
        # So, regardless of where you are in the channel,
        # interpolate all proptable values based on your temp
        if T_m > T_sat:
            T_m = T_sat
        props_Tmz_in = interpPropTable('T (C)', T_m, proptable)

        # calculate h_iplus1, for the first time (this will only happen once)
        if h_iplus1 == False:
            h_i_in = props_Tmz_in['h_l (J/kg)']
            # using .values[0] keeps the index of the value out of our new variable.
        # if h_iplus1 has an actual value (isn't that ~nonsense~ False from earlier)
        # your old h_iplus1 is now your new h_i...
        elif not h_iplus1 == False:
            h_i_in = h_iplus1

        #calculate new h_i+1 based on last CV's h_i and q'_{1+1/2} equation
        h_iplus1 = h_i_in + qp_iPlusHalf(z)*Delta_z/inp.m_dot

        # now interpolate new proptable based on new enthalpy (h_iplus1)
        props_Tmz_out = interpPropTable('h_l (J/kg)', h_iplus1, proptable)

        # take T_m from this new proptable and define it as it's own variable.
        # this is needed to:
        #     (1) export to the final dataframe, 
        #     (2) repeat this calculation up the rest of the channel, and 
        #     (3) compute all other temps.
        T_m = props_Tmz_out['T (C)']
        # UNLESS YOU'RE ABOVE T_sat, THEN USE T_sat PROPERTIES REGARDLESS!!!!!
        if T_m > T_sat:
            T_m = T_sat
        props_Tmz_out = interpPropTable('T (C)', T_m, proptable)
        
        #define properties from props_=_Tmz_out that we need.
        props_Tmz_out_0 = props_Tmz_out.copy()

        rho_l = props_Tmz_out_0['rho_l (kg/m^3)']
        rho_out = rho_l
        rho_g = props_Tmz_out_0['rho_g (kg/m^3)']

        vol_l = props_Tmz_out_0['vol_l (m^3/kg)']
        vol_g = props_Tmz_out_0['vol_g (m^3/kg)']    
        vol_lg = vol_g - vol_l      

        mu_l = props_Tmz_out_0['mu_l (kg/m-s)']
        mu_g = props_Tmz_out_0['mu_g (kg/m-s)']

        h_l = props_Tmz_out_0['h_l (J/kg)']
        h_g = props_Tmz_out_0['h_g (J/kg)']
        h_lg = h_g - h_l

        k_l = props_Tmz_out_0['k_l (W/m-K)']
        Pr = props_Tmz_out_0['Pr_l (arb. unit)']
        # END:   COMPUTE ENERGY BALANCE IN COOLANT Tm(z) #####################



        # BEGIN: COMPUTE X WITH AND WITHOUT SCB ##############################
        # create boiling flags
        Boiling_out_flag = False
        SCB_out_flag = False
        # define the x of the last CV,
        # if this is the first CV, x was previously initialized to 0,
        # so then x_prior == x == 0
        x_prior = x
        # calc x_e from energy balance
        x_e = (h_iplus1 - h_l_sat)/( h_g_sat - h_l_sat)

        # if SCB model is OFF:
        if inp.SCB_flag == False:
            # subcooled: both = 0
            if x_e <= 0:
                x_e = 0
                x = x_e
                SCB_out_flag = False
            # boiling: x = x_e(energy balance)
            else:
                x = x_e
                SCB_out_flag = False

        # if SCB model is ON:
        else:
            if first_zD_flag == True:
                # model is on but T_co is still below T_sat
                x_e = x_e
                x = 0
            else:
                #model is on, and in the last CV, T_co > T_sat    
                x_e = x_e
                # calculate x
                x = x_e - x_e_zD*np.exp(x_e/x_e_zD - 1)           

        # alert user if boiling (this flag is sent to output file)
        if x > 0:
            Boiling_out_flag = True
        # alert user if SCB is occurring (this flag is sent to output file)
        if x > 0 and x_e <= 0:
            SCB_out_flag = True 
        # END:   COMPUTE X WITH AND WITHOUT SCB ##############################

        

        # BEGIN: COMPUTE ONE-PHASE CONVECTIVE HEAT TRANSFER Tco(z) ###########
        #calculate htc using Weisman model:
        # Reynolds number (arb. unit)
        Re = D_e*G/mu_l 
        # Nusselt number, circular tube (arb. unit)
        Nu_ct = 0.023*Re**(0.8)*Pr**(0.333)
        # calculate $\psi$ using Weisman model, a geometry adjustment (arb. unit)
        psi = 1.826*inp.Pitch/inp.D - 1.0430
        Nu = Nu_ct * psi # Nusselt number (arb. unit)
        # calculate Weisman heat transfer coefficient. (W/m^2-K)
        htc_lo = Nu*k_l/D_e

        # BEGIN: COMPUTE TWO-PHASE HEAT TRANSFER T_co(z) #####################
        #calculate q''
        qpp = qp_iPlus1(z)/(2*np.pi*R_co)
        # determine which htc to use based on quality
        if x > 0:
        # if boiling:
            # calculate Xtt, use x and not x_e
            X_tt = ((1-x)/x)**0.9 * (rho_g/rho_l)**0.5 * (mu_l/mu_g)**0.1
            # calculate htc using Schrock and Grossman
            htc_2phase = htc_lo * (7400 * (qpp/G/h_lg) + 1.11*X_tt**-0.66)
            # use this value since we have boiling
            htc_use = htc_2phase
        else:
            # if there's no boiling, use htc_lo
            htc_use = htc_lo
            # calculate T_co using whichever htc is applicable 
        T_co = T_m + qp_iPlus1(z)/(2*np.pi*R_co*htc_use)
        # END:   COMPUTE ONE-PHASE CONVECTIVE HEAT TRANSFER Tco(z) ###########

        # if Tco < sat, SCB is not currently happening
        if T_co < T_sat:
            # if SCB hasn't happened at all, x_e = 0
            if first_zD_flag == True:
                x_e = 0
            # if SCB has happened before, x_e = as usual
            else:
                x_e = x_e
        # if we currently have SCB
        else:
            # if this is the first CV it occurs in
            if first_zD_flag == True:
                # capture this z value
                z_D = z
                # capture this CV's x_e for x calculations in SCB
                x_e_zD = x_e
                # change flag to indicate that any future SCB
                # is NOT the first time it's happening
                first_zD_flag = False
                # output x_e = 0 for this cell, after this, output
                # the actual x_e since SCB has begin
                x_e = 0
            # if you already found z_D, don't change anything
            else:
                z_D = z_D
                x_e_zD = x_e_zD
                x_e = x_e
                first_zD_flag = False
        # END:   COMPUTE TWO-PHASE HEAT TRANSFER T_co(z) #####################



        # BEGIN: COMPTUE CONDUCTION THROUGH CLADDING Tci(z) ##################
        # calculate Tci
        # np.log is ln(), np.log10 is log_base10()
        T_ci = T_co + qp_iPlus1(z)/(2*np.pi*inp.k_cl) * np.log(R_co/R_ci)
        # END:   COMPUTE GAP CONDUCTANCE Tci(z) ##############################



        # BEGIN: COMPUTE GAP CONDUCTANCE T_fo(z) #############################
        # initial heat transfer coefficient guess
        htc_g_guess = 5000 # (W m^-2 K^-1)
        # Temperature, fuel outer
        T_fo_old = 0 #  (C)
        # compute T_fo_new 
        T_fo_new = T_ci + qp_iPlus1(z)/(2*np.pi*R_g*htc_g_guess)
        convergence_flag_Tfo = False

        # while convergence is not achieved 
        while convergence_flag_Tfo == False:
            # Compute k_gas
            # Do this based on avg T in gap in Kelvin using 273.15.
            A = 15.8 # for helium
            T_avg_K = ((T_fo_new + 273.15) + (T_ci + 273.15))/2
            k_gas = A*10**(-4) * T_avg_K**(0.79)
            # compute new htc_g given emmissivities = 1
            htc_g_new = k_gas/delta_eff + inp.SB_constant * ((T_fo_new + 273.15)**4 - (T_ci + 273.15)**4)/((T_fo_new + 273.15) - (T_ci + 273.15))
            # compute T_fo
            T_fo_new = T_ci + qp_iPlus1(z)/(2*np.pi*R_g*htc_g_new)
            # compare old and new T_fo's
            if abs(T_fo_new - T_fo_old) <= inp.convergence_crit_gap:
                # if convergence is achieved, define T_fo and exit the loop.
                convergence_flag_Tfo = True
                T_fo = T_fo_new
            else:
                # if convergence isn't achieved, set _old to equal _new and start over
                convergence_flag_Tfo = False
                T_fo_old = T_fo_new
        # END:   COMPUTE GAP CONDUCTANCE T_fo(z) #############################



        # BEGIN: COMPTUE GAP CONDUCTANCE T_max(z) ############################
        k_bar_guess = 3 # (W m^-1 K^-1)
        #calculate T_max
        T_max_new = T_fo + qp_iPlus1(z)/(4*np.pi*k_bar_guess) 
        # compute conductivity integral at Tmax
        kdT_max = 3824*np.log(402.4 + T_fo) + 6.1256E-11/4 * ((T_fo+273)**4) + qp_iPlus1(z)/(4*np.pi)
        convergence_flag_Tmax = False

        # while convergence is not achieved
        while convergence_flag_Tmax == False:
            kdT_check = 3824*np.log(402.4 + T_max_new) + 6.1256E-11/4 * ((T_max_new+273)**4)
            # demote current T_max_new to "old"
            T_max_old = T_max_new
            # compute revised T_max
            T_max_new = T_max_old + (kdT_max - kdT_check)/100
            #check for convergence
            if abs(kdT_check - kdT_max) <= inp.convergence_crit_fuel:
                convergence_flag_Tmax = True
                T_max = T_max_new
            else:
                convergence_flag_Tmax = False
        # END:   COMPUTE GAP CONDUCTANCE T_max(z) ##########################



        # BEGIN: COMPUTE CHF AND CHFR ########################################
        # Use Bowring correlation with same \psi as Weisman correlation
        P_r = 0.145*(inp.P_nom/1E6) #P must be in MPa
        n = 2 - 0.5*P_r

        # calcualte F_#'s based on P_r
        P_r_exp = np.exp((1-P_r))
        if x > 0:
            if P_r > 1:
                F_1 = P_r**-0.368 * P_r_exp**0.648
                F_3 = P_r**0.219
                F_4 = F_3 * P_r**1.649
                F_2 = F_1 / (P_r**-0.448 * P_r_exp**0.245)
            elif P_r < 1:
                P_r_exp_2089 = P_r_exp**20.89
                P_r_exp_2444 = P_r_exp**2.444
                P_r_exp_16658 = P_r_exp**16.658
                F_1 = (P_r**18.942 * P_r_exp_2089 + 0.917) / 1.917
                F_3 = (P_r**17.023 * P_r_exp_16658 + 0.667) / 1.667
                F_4 = F_3 * P_r**1.649
                F_2 = 1.309 * F_1 / (P_r**1.316 * P_r_exp_2444 + 0.309)
            else:
                F_1 = F_2 = F_3 = F_4 = 1

            A = (2.317*h_lg*D_e*G/4*F_1)/(1 + 0.0143*F_2*G*np.sqrt(D_e))
            B = G*D_e/4
            C = (0.077*F_3*D_e*G)/(1 + 0.347*F_4*(G/1356)**n)
            qpp_CHF = (A - B*h_lg*x)/C * psi # MULTIPLY BY \psi FORM Weisman correlation

            #calculate CHF ratio relative to your current q''
            CHFR = qpp_CHF/qpp

            # if you're over the CHF limit, change the flag, otherwise, don't
            if CHFR < inp.CHFR_crit_limit:
                CHFR_crit_flag = True
            else:
                CHFR_crit_flag = False
        # else: if you're not boiling, you don't have CHF...
        else:
            CHFR = 0
        # END:   COMPUTE CHF AND CHFR ########################################



        # BEGIN: COMPUTE PRESSURE LOSS FOR ONE-PHASE COOLANT #################
        # must use Cheng and Todreas correlation.
        # FOR ONE PHASE: NEGLECT ACCELERATION.
        # if there's no boiling, use one-phase model
        if x <= 0:
            f_interior = Re**(-0.18) * (0.1339 + 0.09059*(inp.Pitch/inp.D - 1) - 0.09926*(inp.Pitch/inp.D - 1)**2)
            # rho*v^2/2 = G^2/(2*rho); inplement below.
            dP_dz_fric = f_interior/D_e*G**2/(2*rho_l)
            # multiple by CV length to get total drop (dP/dz * z = dP)
            DeltaP_fric = dP_dz_fric * Delta_z
            #now, gravitational pressure loss is needed.
            DeltaP_grav = inp.g*(rho_out*Delta_z)
            #add the two together. this is for one-phase only, so no acceleration
            DeltaP_total = DeltaP_fric + DeltaP_grav
        # END: COMPUTE PRESSURE LOSS FOR ONE-PHASE COOLANT #################



        # BEGIN: COMPUTE PRESSURE LOSS FOR TWO-PHASE COOLANT  ##############
        # if there IS boiling, two-phase pressure loss calculated with HEM
        else:
            rho_m = 1/(x/rho_g + (1-x)/rho_l)
            # use same f as single-phase. f_interior = f_lo = f_TP
            f_TP = Re**(-0.18) * (0.1339 + 0.09059*(inp.Pitch/inp.D - 1) - 0.09926*(inp.Pitch/inp.D - 1)**2) 
            dx_dz = (x - x_prior)/Delta_z
            # use HEM, denominator = 1 due to imcompressiblity assumption
            dp_dz_HEM = f_TP/D_e*G**2/(2*rho_m) + G**2*dx_dz*vol_lg + rho_m*inp.g
            DeltaP_total = dp_dz_HEM * Delta_z
        # END:   COMPUTE PRESSURE LOSS FOR TWO-PHASE COOLANT #################

        # put data into the previously defined variable
        data.append([z,
                    T_m,
                    T_co,
                    T_ci,
                    T_fo,
                    T_max,
                    x,
                    x_e,
                    Boiling_out_flag,
                    SCB_out_flag,
                    CHFR,
                    CHFR_crit_flag,
                    DeltaP_total,
                    Re
                    ])
    # END:   OUTER LOOP OF Z VECTOR ##########################################



    # BEGIN: ORGANIZE AND EXPORT DATA ########################################
    # Convert data to a list of dictionaries for easier use
    data_dicts = []
    column_names = ['z', 'T_m (C)', 'T_co (C)', 'T_ci (C)', 'T_fo (C)', 'T_max (C)', 'x', 'x_e', 'Boiling_out_flag', 'SCB_out_flag', 'CHFR', 'CHFR_crit_flag', 'DeltaP_thisCell (Pa)', 'Re']
    for row in data:
        data_dicts.append(dict(zip(column_names, row)))

    # Calculate the sum of DeltaP_thisCell (Pa)
    DeltaP_sum = sum(d['DeltaP_thisCell (Pa)'] for d in data_dicts)

    # Add the DeltaP_sum (Pa) to each dictionary in the list
    for d in data_dicts:
        d['DeltaP_sum (Pa)'] = DeltaP_sum

    timestr = time.strftime('%Y%m%d_%H%M%S')
    dir_output = f'.\\outputs\\SCAcode_output_{infile_name}_{timestr}\\'
    filename = f'SCAcode_data_{infile_name}_{timestr}.csv'
    if not os.path.exists(dir_output):
        os.makedirs(dir_output)
    fp_data = os.path.join(dir_output, filename)

    # Write data to a CSV file
    with open(fp_data, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(data_dicts[0].keys())) # use the keys of the first dictionary as fieldnames
        writer.writeheader()
        for d in data_dicts:
            writer.writerow(d)

    time_end = time.time()
    time_elapsed = time_end - time_start

    print(f"""
{os.path.basename(__file__)} successfully ran.
Elapsed run time: {time_elapsed:.2f} s.
Results were exported to {fp_data}.
    """)
    return fp_data, dir_output
    # END:   ORGANIZE AND EXPORT DATA ########################################
# END:   SOLVER FUNCTION #####################################################