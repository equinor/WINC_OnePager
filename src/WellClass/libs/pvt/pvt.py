from typing import Union, Dict, Tuple, Callable
from pathlib import Path
import numpy as np
import pandas as pd
import json
from scipy.interpolate import RectBivariateSpline
from scipy.integrate import solve_ivp

import scipy.constants as const


def file_exists(files_path: list):
    '''
    Receives a list of file paths and verify if the files in such
    path exists. Throws an exception in case the file is not found
    '''
    for file in files_path:
        if not file.exists():
            raise FileNotFoundError(f"Required PVT file not found: {file}")


def get_mixture_info(pvt_path: Union[str, Path]):
    with open(f"{pvt_path}/metadata.json", "r") as file:
        mixture_info = json.load(file)
    return mixture_info


def corr_rhobrine_LaliberteCopper(salinity: float, temperature: np.ndarray, pressure: np.ndarray, rho_h2o: np.ndarray) -> np.ndarray:
    '''
    Reads the vectors for pressure and temperature and the matrix for rho
    Note that the values for temperature and rho must be aligned with values in rho
    Note also for rho: One temperature for each column
                        One pressure for each row
    It applies a salinity correction for a concentration of X gNaCl in 100 g H2O
    based on the work done by 
    LALIBERTÉ, M. & COOPER, W. E. 2004. Model for Calculating the Density of 
    Aqueous Electrolyte Solutions. Journal of Chemical & Engineering Data, 49, 
    1141-1151. https://doi.org/10.1021/je0498659
    https://www.calsep.com/13-density-of-brine/
    '''


    ## Laliberté and Cooper model for NaCl solutions
    # Laliberté and Cooper model: constants for NaCl
    c0 = -0.00433
    c1 =  0.06471
    c2 = 1.0166
    c3 = 0.014624
    c4 = 3315.6

    # NaCl concentration
    w = salinity / 100

    # Laliberté and Cooper model: Apparent density 
    rho_app = (c0*w + c1)*np.exp(0.000001 * (temperature + c4)**2) / (w + c2 + c3 * temperature)

    # Laliberté and Cooper model: Brine density
    rho_brine =  1 / (((1-w)/rho_h2o) + (w/rho_app))

    return rho_brine


def load_pvt_data(pvt_root_path: Union[str, Path], fluid_type: str = None, load_fluid: bool = True) -> Dict[str, Dict[str, np.ndarray]]:
    # Convert to Path object if not already
    pvt_root_path = Path(pvt_root_path)

    # Define file names for shared PVT data files
    temperature_file = pvt_root_path / "temperature.txt"
    pressure_file = pvt_root_path / "pressure.txt"
    
    # Define file names for brine data files
    brine_path = pvt_root_path / "water"
    brine_metadata_file = brine_path / "metadata.json"
    brine_rho_file = brine_path / "rho.txt"

    # Check that shared files exist before proceeding
    file_exists([temperature_file, pressure_file, brine_metadata_file, brine_rho_file])

    # Load shared data from files
    temperature_data = np.loadtxt(temperature_file)
    pressure_data = np.loadtxt(pressure_file)

    # Load density data
    rho_brine = np.loadtxt(brine_rho_file, delimiter=',')  # Assuming comma separation

    # Load metadata for brine
    brine_metadata = get_mixture_info(brine_path)

    # Create a dictionary for brine PVT data
    pvt_data = {
        "temperature": temperature_data,
        "pressure": pressure_data,
        "brine": {
            "rho": rho_brine,
            "metadata": brine_metadata,
        }
    }

    # If reservoir fluid data is requested, load it
    if load_fluid and fluid_type is not None:

        # Load fluid-specific data
        fluid_path = pvt_root_path / fluid_type
        fluid_metadata_file = fluid_path / "metadata.json"
        fluid_rho_file = fluid_path / "rho.txt"

        # Check that fluid-specific files exist before proceeding
        file_exists([fluid_metadata_file, fluid_rho_file])

        # Load density data for the specified fluid
        rho_fluid = np.loadtxt(fluid_rho_file, delimiter=',')  # Assuming comma separation

        # Load metadata for the specified fluid
        fluid_metadata = get_mixture_info(fluid_path)

        # Add fluid data to the PVT data dictionary
        pvt_data[fluid_type] = {
            "rho": rho_fluid,
            "metadata": fluid_metadata 
        }

    return pvt_data


def compute_hydrostatic_pressure(depth_array: np.ndarray, temperature_array: np.ndarray, 
                                 pvt_data: Dict[str, np.ndarray], salinity: float = 3.5) -> np.ndarray:
    """
    Compute the hydrostatic pressure profile based on depth and temperature arrays and PVT data for water.

    Args:
        depth_array (np.ndarray): Array of depth values.
        temperature_array (np.ndarray): Array of temperature values corresponding to the depth values.
        pvt_data (Dict[str, np.ndarray]): Dictionary containing PVT data, including pressure and water density.

    Returns:
        np.ndarray: Array of hydrostatic pressure values corresponding to the depth values.
    """
    # Extract pressure and water density data from PVT data
    p_vec = pvt_data["pressure"]
    t_vec = pvt_data["temperature"]
    rho_h2o_vec = pvt_data['brine']["rho"]

    t_grid, p_grid = np.meshgrid(p_vec, t_vec)

    # Correct water density for salinity using Laliberté and Cooper model
    rho_brine = corr_rhobrine_LaliberteCopper(salinity, t_grid, p_grid, rho_h2o_vec)

    # Create an interpolator for water density based on PVT data
    get_rho_brine = RectBivariateSpline(p_vec, pvt_data["temperature"], rho_brine)

    # Define the ODE system for hydrostatic pressure integration
    def odesys(z, P):
        T = np.interp(z, depth_array, temperature_array)  # Interpolate temperature at depth z
        rho = get_rho_brine(P, T)  # Get water density at pressure P and temperature T
        dPdz = rho * const.g / const.bar  # Pressure gradient
        return dPdz

    # Initial conditions: atmospheric pressure at the surface
    P_0 = const.atm / const.bar

    # Solve ODEs from the surface to the final depth
    solution = solve_ivp(odesys, [depth_array[0], depth_array[-1]], [P_0], t_eval=depth_array, vectorized=True, method='Radau')

    # Return the hydrostatic pressure profile
    return solution.y[0]


# Define _Pdz_odesys outside with necessary arguments
def _Pdz_odesys(z: float, P: float, 
                depth_array: np.ndarray, 
                temperature_array: np.ndarray, 
                rho_interpolator: RectBivariateSpline, 
                fluid_key: str) -> Tuple[float]:
    
    T = np.interp(z, depth_array, temperature_array)
    rho = get_rho_from_pvt_data(P, T, rho_interpolator)
    dPdz = rho * const.g / const.bar
    return dPdz


def _integrate_pressure(init_curves: pd.DataFrame, 
                        reference_depth: float, 
                        reference_pressure: float,
                        pvt_data: dict, 
                        fluid_key: str,
                        interpolator: RectBivariateSpline,
                        top_limit: float = None, bottom_limit:float = None) -> np.ndarray:
    '''
    Simple integration to find pressure
    Starting point: reference_depth at reference pressure.  Reference temperature is found in init_curves
                    Then iterate upwards (up) or downwards (down) to subtract or add pressure.
                    Recalculates rho at each step.
    '''
    depth_array = init_curves['depth'].values

    ## Initialization
    #New columns needed in DataFrame
    if fluid_key == 'brine':
        fluid_label = fluid_key
    else:
        fluid_label = 'fluid'

    colname_p = fluid_label+'_pressure'


    if colname_p not in init_curves.columns:
        init_curves[colname_p] = np.nan

    # Check and update limits
    if top_limit is None:
        top_limit = float(depth_array.min())
    
    if bottom_limit is None:
        bottom_limit = float(depth_array.max())

    # Check if reference_depth is within limits
    if reference_depth < top_limit or reference_depth > bottom_limit:
        raise ValueError('reference_depth is out of range')
    
    # Check if the reference_depth value is equal to bottom AND top limit
    if reference_depth == top_limit and reference_depth == bottom_limit:
        raise ValueError('reference_depth cannot be equal to bottom and top limit')

    # Check and integrate upwards if needed
    if reference_depth >= top_limit:
        ivp_data = build_ivp_data(reference_depth, top_limit, reference_pressure, init_curves)
        ivp_solution = solve_ivp_with_data(ivp_data, interpolator, fluid_key)
        init_curves.loc[ivp_data["t_eval"].index, colname_p] = ivp_solution.y[0]

    # Check and integrate downwards if needed
    if reference_depth <= bottom_limit:
        ivp_data = build_ivp_data(reference_depth, bottom_limit, reference_pressure, init_curves)
        ivp_solution = solve_ivp_with_data(ivp_data, interpolator, fluid_key)
        init_curves.loc[ivp_data["t_eval"].index, colname_p] = ivp_solution.y[0]

    if reference_depth in init_curves['depth'].values:
        init_curves.loc[init_curves['depth'] == reference_depth, colname_p] = reference_pressure

    return init_curves[colname_p].values.astype(float)


def build_ivp_data(interval_start: float,
                   interval_end: float,
                   initial_state: float,
                   init_curves: pd.DataFrame) -> Dict:
    """
    Builds a dictionary containing all the necessary values to solve the IVP

    Args:
        interval_start (float): The start reference depth
        interval_end   (float): The end reference depth
        initial_state  (float): Reference pressure value
        init_curves    (pd.Dataframe): Data frame containing 
                                        the data values

    Returns:
        Dict: A dictionary containing the necessary VALUES to solve the IVP. Considering that the solve_ivp
            method has the following signature: 
                solve_ivp(fun, t_span, y0, method='RK45', t_eval=None, dense_output=False, events=None, 
                    vectorized=False, args=None, **options)
        
        The following can be used:
            t_span: [ivp_data["interval_start"], ivp_data["interval_end"]]
            y0: [ivp_data["initial_state"]], 
            t_eval = ivp_data["t_eval"].values,
            args = (ivp_data["depth_array"].values, 
                    ivp_data["temperature_array"].values
                    /*other arguments not related to ivp_data*/),

        This method DOES NOT account for:
         1. The function to solve the IVP 
         2. Default parameters of the 'solve_ivp' method (e.g.: method, dense_output)
         3. Arguments that have no need to be calculated or extracted from the input data frame
    """
    
    query = None
    if (interval_start >= interval_end):
        query = init_curves[init_curves.depth < interval_start]
        query = query.sort_values('depth', ascending = False)
    elif (interval_start <= interval_end):
        query = init_curves[init_curves.depth > interval_start]

    ivp_data = {
        "interval_start": interval_start,
        "interval_end": interval_end,
        "initial_state": initial_state,
        "t_eval": query['depth'],
        "depth_array": init_curves['depth'],
        "temperature_array": init_curves['temperature']
    }
    return ivp_data


def solve_ivp_with_data(ivp_data: dict,
                        interpolator: RectBivariateSpline,
                        fluid_key: str):
    
    ivp_solution = solve_ivp(fun = _Pdz_odesys,
                             t_span = [ivp_data["interval_start"], ivp_data["interval_end"]], 
                             y0 = [ivp_data["initial_state"]], 
                             t_eval = ivp_data["t_eval"].values,
                             dense_output = True, 
                             args = (ivp_data["depth_array"].values,
                                      ivp_data["temperature_array"].values, 
                                      interpolator, 
                                      fluid_key),
                             method = 'Radau')
    return ivp_solution


def get_rho_from_pvt_data(pressure: float, temperature: float, rho_interpolator: RectBivariateSpline) -> float:
    # Interpolate the density at the given pressure and temperature
    rho = rho_interpolator(pressure, temperature)[0, 0]
    return rho