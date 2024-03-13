
import os
import numpy as np
import pandas as pd

import scipy
import scipy.constants
from scipy.interpolate import RectBivariateSpline

'''Some global parameters'''
G       = scipy.constants.g   #9.81 m/s2 gravity acceleration

def get_pvt(pvt_path: str) -> tuple:
    '''Reads the vectors for pressure and temperature and the matrix for rho
       Note that the values for temperature and rho must be aligned with values in rho
       Note also for rho: One temperature for each column
                          One pressure for each row
       It applies a salinity correction for a concentration of 3.5 gNaCl in 100 g H2O
       based on the work done by 
       LALIBERTÉ, M. & COOPER, W. E. 2004. Model for Calculating the Density of 
       Aqueous Electrolyte Solutions. Journal of Chemical & Engineering Data, 49, 
       1141-1151.
       https://www.calsep.com/13-density-of-brine/
    '''
    fn_temp    = os.path.join(pvt_path, "temperature.txt")
    fn_pres    = os.path.join(pvt_path, "pressure.txt")
    fn_rho_co2 = os.path.join(pvt_path, "rho_co2.txt")
    fn_rho_h2o = os.path.join(pvt_path, "rho_h2o.txt")

    t   = np.loadtxt(fn_temp)
    p   = np.loadtxt(fn_pres)
    rho_co2 = np.loadtxt(fn_rho_co2,delimiter=',')
    rho_h2o = np.loadtxt(fn_rho_h2o,delimiter=',')

    #compute 2d matrices for pressure and temperature
    t_grid, p_grid = np.meshgrid(t, p)

    # Laliberté and Cooper model for NaCl solutions
    # Laliberté and Cooper model: constants for NaCl
    c0 = -0.00433
    c1 =  0.06471
    c2 = 1.0166
    c3 = 0.014624
    c4 = 3315.6

    # NaCl concentration
    w = 3.5 / 100

    # Laliberté and Cooper model: Apparent density 
    rho_app = (c0*w + c1)*np.exp(0.000001 * (t_grid + c4)**2) / (w + c2 + c3 * t_grid)

    # Laliberté and Cooper model: Brine density
    rho_brine =  1 / (((1-w)/rho_h2o) + (w/rho_app))

    return t, p, rho_co2, rho_brine

def get_hydrostatic_P(well_header: dict, *, dz=1, pvt_path: str) -> pd.DataFrame:
    '''Simple integration to get the hydrostatic pressure at a given depth
       Does also calculates the depth column, temperatur vs depth and water density (RHOH2O) vs depth (hydrostatic)
    '''
    t_vec, p_vec, rho_co2_vec, rho_h2o_vec = get_pvt(pvt_path)

    #Make interpolators for the imported tables
    get_rho_h2o = RectBivariateSpline(p_vec, t_vec, rho_h2o_vec)

    #Make the depth-vector from msl and downwards
    td_msl = well_header['well_td_rkb']-well_header['well_rkb']
    z_vec  = np.arange(0, int(td_msl)+500, dz)

    #Create dataframe for storing pressures and temperatures. hs_p_df -> HydroStatic_Pressure_DataFrame
    hs_p_df = pd.DataFrame(data=z_vec, columns = ['depth_msl'])

    #Compute temperature. Constant in water column and as a function of input geothermal gradient
    hs_p_df['temp'] = well_header['sf_temp'] + (hs_p_df['depth_msl']-well_header['sf_depth_msl'])*(well_header['geo_tgrad']/1000)
    hs_p_df.loc[hs_p_df['depth_msl']<well_header['sf_depth_msl'], 'temp'] = well_header['sf_temp']


    ##Integrate hydrostatic pressure
    #Pressure (atm), depth and temperature at msl
    p0 = scipy.constants.atm/scipy.constants.bar  #1.01325 bar Pressure at MSL
    z0 = 0
    t0 = np.interp(z0, hs_p_df['depth_msl'], hs_p_df['temp'])

    #Start to integrate downwards from msl. Assign first entry at zero depth (msl).
    rho_vec = [get_rho_h2o(p0, t0)[0,0]]
    hs_p_df['hs_p'] = p0 

    p = p0
    for idx, t in hs_p_df['temp'][1:].items():
        rho = get_rho_h2o(p,t)[0,0]
        p += (rho*G*dz)/scipy.constants.bar    #dp = rho*g*h/1e-5  the latter to go from Pascal to atm
        rho_vec.append(rho)
        hs_p_df.loc[idx, 'hs_p'] = p
    hs_p_df['RHOH2O'] = rho_vec

    return hs_p_df

