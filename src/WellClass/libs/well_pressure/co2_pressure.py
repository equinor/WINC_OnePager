
import numpy as np
import pandas as pd

from scipy.interpolate import RectBivariateSpline
import scipy.constants as const
from typing import Union, Tuple
from scipy.integrate import solve_ivp

from ..pvt.pvt import get_hydrostatic_P, get_pvt
from ..utils.compute_intersection import compute_intersection

'''Some global parameters'''
G       = const.g   #9.81 m/s2 gravity acceleration
BAR2PA  = const.bar #10**5 Going from bars to Pascal: 1 bar = 10**5 Pascal
SHMIN_FAC = 0.1695

'''Global names. Should be replaced by a global class containing these names? TODO'''
SF_DEPTH_NAME = 'sf_depth_msl'
DEPTH_NAME    = 'depth_msl'
SHMIN_NAME    = 'Shmin'
TEMP_NAME     = 'temp'
RHO_NAME      = 'rho'
MAX_PRESSURE_NAME = 'max_pressure'


def _get_shmin(well_header: dict, pt_df_in: pd.DataFrame) -> pd.DataFrame:
    '''Compute Shmin based on an empirical formula
    '''
    pt_df = pt_df_in.copy()
    sf_depth = well_header['sf_depth_msl']                                 #Depth to sea floor
    sf_pressure = np.interp(sf_depth, pt_df['depth_msl'], pt_df['hs_p'])   #Pressure at sea floor
    pt_df[SHMIN_NAME] = pt_df['hs_p']                                         #Initialize with hydrostatic pressure
    shmin_query = pt_df.query('depth_msl>=@sf_depth')                      #Take only values below sea floor
    pt_df.loc[pt_df['depth_msl']>=sf_depth, SHMIN_NAME] = sf_pressure + (shmin_query['depth_msl']-sf_depth)*SHMIN_FAC

    return pt_df


def _get_rho_in_pressure_column(pt_df_in: pd.DataFrame, colname_p: str, colname_rho: str, get_rho: callable) -> pd.DataFrame:
    '''
    Calculate the water density in the CO2 column with its calculated pressure and temperature.
    
    Parameters:
    - pt_df_in (pd.DataFrame): Input DataFrame containing temperature and pressure data.
    - colname_p (str): Column name in the DataFrame that contains pressure values.
    - colname_rho (str): Column name to be added to the DataFrame for storing calculated density values.
    - get_rho (callable): Function that takes pressure and temperature as inputs and returns the density.
    
    Returns:
    - pd.DataFrame: DataFrame with an additional column for the calculated densities.
    '''
    # Create a copy of the input DataFrame to avoid modifying the original data
    pt_df = pt_df_in.copy()
    
    # Extract temperature and pressure columns from the DataFrame
    temps = pt_df['temp']
    pressures = pt_df[colname_p]
    
    # Initialize an empty list to store the calculated densities
    rhos = []
    
    # Loop through each pair of temperature and pressure values
    for t, p in zip(temps, pressures):
        # Call the get_rho function to calculate the density and append the result to the list
        rhos.append(get_rho(p, t)[0, 0])
    
    # Add the calculated densities as a new column in the DataFrame
    pt_df[colname_rho] = rhos[:]
    
    return pt_df

def _get_max_pressure(pt_df_in: pd.DataFrame, max_pressure_pos: Union[dict, list, float, int], get_rho: callable) -> pd.DataFrame:
    '''
    Calculates downwards maximum pressure at a depth deeper than the input depth.
    I.e. the pressure at deeper locations such that witha CO2-column the pressure will be equal to Shmin at the "start/input"-depth
    Input:
      max_pressure_pos is a depth from where max pressure (wrt Shmin) is calculated. It can be a dict (my_well.barriers) a list of numbers or scalars.
      If my_well.barriers is given then it is calculated from the base of each barrier
    '''

    pt_df = pt_df_in.copy()
    if isinstance(max_pressure_pos, dict):   #Then max_pressure_pos is the same as barriers - and max pressure is calcualted for each barrier
        print(f'max_pressure_pos is a dictionary of barrriers')
        for idx, key in max_pressure_pos['barrier_name'].items():
            barr_depth = max_pressure_pos['bottom_msl'][idx]
            colname_p = f"{MAX_PRESSURE_NAME}_{key}"
            print(f"Calculating max pressure below barrier {key} from depth {barr_depth}")
            p0 = np.interp(barr_depth, pt_df['depth_msl'], pt_df[SHMIN_NAME])
            
            pt_df = _integrate_pressure(pt_df, get_rho, barr_depth, p0, 'down', colname_p)
    elif isinstance(max_pressure_pos, (list, float, int)):
        print(f'max_pressure_pos is a value')
        if isinstance(max_pressure_pos, (float, int)): #Make it a list with one element
            max_pressure_pos = [max_pressure_pos]
        for depth in max_pressure_pos:

            colname_p = f"{MAX_PRESSURE_NAME}_at_{int(depth)}" 
            print(f"Calculating max pressure from depth {depth}")            
            p0 = np.interp(float(depth), pt_df['depth_msl'], pt_df[SHMIN_NAME])
            pt_df = _integrate_pressure(pt_df, get_rho, float(depth), p0, 'down', colname_p)

    return pt_df




def _Pdz_odesys(z: float, y: np.ndarray, well_header: dict, temp_getter: callable, rho_getter: callable)  -> Tuple[float]:
    """
    Calculate the pressure gradient over a depth interval in a well.

    This function computes the derivative of pressure with respect to depth (dPdz) in a fluid column 
    given the depth and current state of pressure. It can be used within numerical methods to solve 
    ordinary differential equations representing the pressure profile along the well.

    Parameters:
    - z (float): The current depth in the well at which we are calculating the pressure gradient.
    - y (np.ndarray): A numpy array containing current state variables, in this context only pressure P.
    - well_header (dict): A dictionary containing metadata or parameters of the well, used by the
                          temperature getter function to compute the current temperature.
    - temp_getter (callable): A function which returns the temperature T at a given depth z 
                              and using well parameters from well_header.
    - rho_getter (callable): A function which returns the density rho of the fluid as a 2D numpy array 
                             at the provided pressure P and temperature T.
    
    Returns:
    - Tuple[float]: A tuple containing a single element, the pressure gradient dPdz (in bar/meter) at depth z.

    """

    P = y[0]
    T = temp_getter(z, well_header)
    rho = rho_getter(P, T)[0, 0]
    dPdz = rho * const.g / const.bar
    return dPdz,


# Compute the temperature given the input gradient
def compute_T(z : Union[float, int] , well_header: dict) -> float:
    """
    Compute the temperature at a given depth in a well using a linear geothermal gradient.
    
    This function calculates the temperature at a specified depth based on the surface
    temperature and the geothermal gradient provided in the well header data. It assumes
    a linear relationship between temperature and depth beyond a reference depth measure.

    Parameters:
    - z (Union[float, int]): The depth at which to compute the temperature, which can be
                             either a float or an integer value.
    - well_header (dict): A dictionary containing specific parameters of the well that
                          are used to calculate the temperature. The expected keys of
                          the dictionary are:
                          'sf_temp' - Surface temperature at seafloor (in degrees Celsius).
                          'sf_depth_msl' - Seafloor depth (mean sea level) in meters where 'sf_temp' is measured.
                          'geo_tgrad' - Geothermal temperature gradient (in degrees Celsius per kilometer).

    Returns:
    - float: The computed temperature in degrees Celsius at the given depth z.
    """

    T = well_header['sf_temp'] + max(0, z - well_header['sf_depth_msl']) * (well_header['geo_tgrad'] / 1000)
    return T


def _integrate_pressure(well_header: dict, pt_df_in: pd.DataFrame, get_rho: callable, reference_depth: float, 
                        reference_pressure: float, colname_p: str, 
                        top_limit: float = None, bottom_limit:float = None) -> pd.DataFrame:
    '''
    Simple integration to find pressure
    Starting point: reference_depth at reference pressure.  Reference temperature is found in pt_df
                    Then iterate upwards (up) or downwards (down) to subtract or add pressure.
                    Recalculates rho at each step.
    '''
    pt_df = pt_df_in.copy()

    print(f'{reference_pressure=}, {reference_depth=}')

    ## Initialization
    #New columns needed in DataFrame
    if colname_p not in pt_df.columns:
        pt_df[colname_p] = np.nan

    colname_rho = colname_p+'_rho'

    if colname_rho not in pt_df.columns:
        pt_df[colname_rho] = np.nan

    # Check and update limits
    if top_limit is None:
        top_limit = float(pt_df.depth_msl.min())
    
    if bottom_limit is None:
        bottom_limit = float(pt_df.depth_msl.max())

    # Check if reference_depth is within limits
    if reference_depth < top_limit or reference_depth > bottom_limit:
        raise ValueError('reference_depth is out of range')

    # Check and integrate upwards if needed
    if reference_depth > top_limit:
    
        query_up = pt_df[pt_df.depth_msl < reference_depth]
        query_up = query_up.sort_values('depth_msl', ascending=False)


        solution_up = solve_ivp(_Pdz_odesys, [reference_depth, top_limit], [reference_pressure], args=(well_header, compute_T, get_rho), 
                         t_eval=query_up['depth_msl'].values, dense_output=True,
                         method = 'Radau')
    
        pt_df.loc[query_up.index, colname_p] = solution_up.y[0]

    # Check and integrate downwards if needed
    if reference_depth < bottom_limit:

        query_down = pt_df[pt_df.depth_msl > reference_depth]
        solution_down = solve_ivp(_Pdz_odesys, [reference_depth, bottom_limit], [reference_pressure], args=(well_header, compute_T, get_rho), 
                         t_eval=query_down['depth_msl'].values, dense_output=True,
                         method = 'Radau')

        pt_df.loc[query_down.index, colname_p] = solution_down.y[0]


    # Vectorized retrieval of densities
    query_total = pt_df.query('depth_msl>=@top_limit and depth_msl<=@bottom_limit')
    P_array = pt_df.loc[query_total.index, colname_p].values
    T_array = pt_df.loc[query_total.index, 'temp'].values

    densities = np.array([get_rho(P, T) for P, T in zip(P_array, T_array)])

    # After integration, access the stored densities in dataframe
    pt_df.loc[query_total.index, colname_rho] = densities.flatten()

    

    return pt_df

def compute_MSAD(p_init: dict, pt_df: pd.DataFrame):

    """
    Calculates MSAD: Minimum Safety Abandonement Depth
    MSAD is the intersection point between the CO2 pressure and Shmin.
    It computes a pressure and depth value for every pressure scenario.
    """

    MSAD = dict()

    shmin = pt_df.Shmin.values
    depth = pt_df.depth_msl.values


    for key, value in p_init.items():
        if key == 'depth_msl':
            print(f"Reference depth: {value}")

        else:

            MSAD[key] = dict()

            co2_p = pt_df[f'{key}_co2'].values
            z_MSAD, p_MSAD = compute_intersection(x = depth, y1 = shmin, y2 = co2_p)

            MSAD[key]['z_MSAD'] = z_MSAD
            MSAD[key]['p_MSAD'] = p_MSAD


    return MSAD


    




