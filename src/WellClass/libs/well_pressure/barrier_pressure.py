
import numpy as np
import pandas as pd

import scipy.constants

'''Some global parameters'''
G       = scipy.constants.g   #9.81 m/s2 gravity acceleration
BAR2PA  = scipy.constants.bar #10**5 Going from bars to Pascal: 1 bar = 10**5 Pascal
REGR_A  = -0.000116           #Intercept term in regression equation for the proxy. Consider as a input
REGR_B  = 0.000002725         #Inclination-term in regression equation for the proxy

def compute_barrier_leakage(barrier_perm: dict, reservoir_P: dict, pressure_CO2: dict, barrier_props: dict) -> pd.DataFrame:
    """ compute leakage from the given barrier
    """

    # Calculates pressure above and below the barrier and densities below the barrier
    barrier_p_rho = _get_barrier_p_and_rho(reservoir_P, pressure_CO2, barrier_props)

    # Estimate CO2 leakage in [tons/day] after a trancient period
    barrier_leakage = _get_barrier_leakage(barrier_perm, barrier_p_rho, barrier_props)

    return barrier_leakage

def _get_barrier_p_and_rho(reservoir_P: dict, pressure_CO2: dict, barrier_props: dict) -> pd.DataFrame:
    ''' Calculates pressure above and below the barrier 
        and densities below the barrier using the assumption 
        that the borehole is filled with water above the barrier 
        and filled with CO2 below the barrier
    '''

    #Get the pressure cases to include
    rp_names = []
    for key in reservoir_P:
        if key.startswith("RP"):
            rp_names.append(key)

    df = pd.DataFrame(columns=["p_h2o_above_barrier", "p_co2_below_barrier", "rho_h2o_below_barrier", "rho_co2_below_barrier"], index=rp_names)

    #
    depth = pressure_CO2["depth_msl"]           #Depth look up table used in the interpolation of pressure and densities in the pressure_CO2 dataframe.
    top   = barrier_props['top']
    bottom= barrier_props['bottom']

    df["p_h2o_above_barrier"] = [np.interp(top,    depth, pressure_CO2["hydrostatic_pressure_h2o"])]*len(rp_names)
    df["p_co2_below_barrier"] = [np.interp(bottom, depth, pressure_CO2[f"{key}_co2"]) for key in rp_names]

    df["rho_h2o_below_barrier"]  = [np.interp(bottom, depth, pressure_CO2[f"{key}_h2o_rho_in_co2_column"]) for key in rp_names]
    df["rho_co2_below_barrier"]  = [np.interp(bottom, depth, pressure_CO2[f"{key}_co2_rho"]) for key in rp_names]

    # barrier_props[barrier_name]['p_and_rho'] = df.copy()
    return df.copy()

def _get_barrier_leakage(barrier_perm: dict, barrier_p_rho: pd.DataFrame, barrier_props: dict) -> pd.DataFrame:
    ''' Returns an estimate of CO2 leakage in [tons/day] after a trancient period.
        It also is based on than the reservoir pressur does not change - hence if the process is completely stationary except for the leakage, 
        i.e. the leaked volumes are << than the reservoir volumes.
        The proxy is based on two steps:
        proxy1: r*r*k/l*(g*drho*l + dp*10**5)    -> To get the most important physical properties determining rate
        proxy2: a + b*proxy1                     -> To calibrate to actual rates estimated by a large number of runs done in pflotran.

        The variables in the proxy regression must come from somewhere.
        Here it is hardcoded in the header, but one could imagine to have several models fit for different circumstances.
        Then the parameteters could be case dependent and an input
    '''

    #Get the permeabilty values to use
    perms = barrier_perm['kv'].values()

    #Get the pressure cases to use (RP1, RP2 etc)
    cases = barrier_p_rho.index

    #Make data-structurs ready.
    df_leakage = pd.DataFrame(columns=perms, index=cases)

    #
    df_p_rho = barrier_p_rho.copy()

    #To make the formula below easier to read
    r      = barrier_props['radius']
    length = float(barrier_props['height'])

    for k in perms:                                                         #Loop the permeability-cases
        for case in cases:                                                  #Loop the pressure cases
            drho = df_p_rho.loc[case, 'rho_h2o_below_barrier'] - df_p_rho.loc[case, 'rho_co2_below_barrier']
            dp   = df_p_rho.loc[case, 'p_co2_below_barrier']   - df_p_rho.loc[case, 'p_h2o_above_barrier']

            prox = (r*r*k/length)*(G*length*drho + dp*BAR2PA)
            df_leakage.loc[case,k] = np.round(max(REGR_A + REGR_B*prox,0),5)
    
    print(df_leakage)

    # barrier_props[barrier_name]['leakage'] = df_leakage.copy()
    return df_leakage.copy()