from typing import Union, Dict, Tuple, Callable
from pathlib import Path
import numpy as np
import pandas as pd
import json
from scipy.interpolate import RectBivariateSpline
from scipy.integrate import solve_ivp

import scipy.constants as const

def load_pvt_data(pvt_root_path: Union[str, Path], fluid_type: str) -> Dict[str, Dict[str, np.ndarray]]:
    # Convert to Path object if not already
    pvt_root_path = Path(pvt_root_path)

    # Define file names for shared PVT data files
    temperature_file = pvt_root_path / "temperature.txt"
    pressure_file = pvt_root_path / "pressure.txt"

    # Check that shared files exist before proceeding
    for file in [temperature_file, pressure_file]:
        if not file.exists():
            raise FileNotFoundError(f"Required PVT file not found: {file}")

    # Load shared data from files
    t = np.loadtxt(temperature_file)
    p = np.loadtxt(pressure_file)

    # Load brine data
    brine_path = pvt_root_path / "water"
    brine_metadata_file = brine_path / "metadata.json"
    brine_rho_file = brine_path / "rho.txt"

    # Load fluid-specific data
    fluid_path = pvt_root_path / fluid_type
    fluid_metadata_file = fluid_path / "metadata.json"
    fluid_rho_file = fluid_path / "rho.txt"

    # Check that fluid-specific files exist before proceeding
    for file in [brine_metadata_file, brine_rho_file, fluid_metadata_file, fluid_rho_file]:
        if not file.exists():
            raise FileNotFoundError(f"Required PVT file not found: {file}")

    # Load density data
    rho_brine = np.loadtxt(brine_rho_file, delimiter=',')  # Assuming comma separation
    rho_fluid = np.loadtxt(fluid_rho_file, delimiter=',')  # Assuming comma separation

    # Load metadata
    with open(brine_metadata_file, 'r') as file:
        brine_metadata = json.load(file)
    with open(fluid_metadata_file, 'r') as file:
        fluid_metadata = json.load(file)

    # Create a dictionary to map each type of data to its corresponding numpy array
    pvt_data = {
        "temperature": t,
        "pressure": p,
        "brine": {
            "rho": rho_brine,
            "metadata": brine_metadata,
        },
        fluid_type: {
            "rho": rho_fluid,
            "metadata": fluid_metadata,
        }
    }

    return pvt_data



def compute_hydrostatic_pressure(depth_array: np.ndarray, temperature_array: np.ndarray, pvt_data: Dict[str, np.ndarray]) -> np.ndarray:
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
    rho_h2o_vec = pvt_data['brine']["rho"]

    # Create an interpolator for water density based on PVT data
    get_rho_h2o = RectBivariateSpline(p_vec, pvt_data["temperature"], rho_h2o_vec)

    # Define the ODE system for hydrostatic pressure integration
    def odesys(z, P):
        T = np.interp(z, depth_array, temperature_array)  # Interpolate temperature at depth z
        rho = get_rho_h2o(P, T)  # Get water density at pressure P and temperature T
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
                get_rho_func: Callable[[float, float, Dict, str], float], 
                pvt_data: Dict, 
                fluid_key: str) -> Tuple[float]:
    T = np.interp(z, depth_array, temperature_array)
    rho = get_rho_func(P, T, pvt_data, fluid_key)
    dPdz = rho * const.g / const.bar
    return dPdz,


def _integrate_pressure(init_curves: pd.DataFrame, 
                        reference_depth: float, 
                        reference_pressure: float,
                        pvt_data: dict, 
                        fluid_key: str,
                        get_rho_func: Callable[[float, float, Dict, str], float],
                        top_limit: float = None, bottom_limit:float = None) -> np.ndarray:
    '''
    Simple integration to find pressure
    Starting point: reference_depth at reference pressure.  Reference temperature is found in init_curves
                    Then iterate upwards (up) or downwards (down) to subtract or add pressure.
                    Recalculates rho at each step.
    '''
    
    
    # ip_arrays = {'depth': depth_array, 'temperature': temperature_array}
    # init_curves = pd.DataFrame(data = ip_arrays)

    depth_array = init_curves['depth'].values
    temperature_array = init_curves['temperature'].values


    # def _Pdz_odesys(z: float, y: np.ndarray)  -> Tuple[float]:
    #     P = y[0]

    #     T = np.interp(z, depth_array, temperature_array)

    #     rho = get_rho_from_pvt_data(P, T, pvt_data, fluid_key)
    #     dPdz = rho * const.g / const.bar
    #     return dPdz,

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

    # Check and integrate upwards if needed
    if reference_depth >= top_limit:
    
        query_up = init_curves[init_curves.depth < reference_depth]
        query_up = query_up.sort_values('depth', ascending=False)

        solution_up = solve_ivp(_Pdz_odesys, [reference_depth, top_limit], [reference_pressure], 
                         t_eval=query_up['depth'].values, dense_output=True,
                         args=(depth_array, temperature_array, get_rho_func, pvt_data, fluid_key),
                         method = 'Radau')
    
        init_curves.loc[query_up.index, colname_p] = solution_up.y[0]

    # Check and integrate downwards if needed
    if reference_depth <= bottom_limit:

        query_down = init_curves[init_curves.depth > reference_depth]
        solution_down = solve_ivp(_Pdz_odesys, [reference_depth, bottom_limit], [reference_pressure], 
                         t_eval=query_down['depth'].values, dense_output=True,
                         args=(depth_array, temperature_array, get_rho_func, pvt_data, fluid_key),
                         method = 'Radau')

        init_curves.loc[query_down.index, colname_p] = solution_down.y[0]

    if reference_depth in init_curves['depth'].values:
        init_curves.loc[init_curves['depth'] == reference_depth, colname_p] = reference_pressure

    return init_curves[colname_p].values.astype(float)


def get_rho_from_pvt_data(pressure: float, temperature: float, pvt_data: Dict[str, Dict[str, np.ndarray]], fluid_key: str = 'brine') -> float:
    # Retrieve the temperature and pressure vectors from the PVT data
    temperature_vector = pvt_data['temperature']
    pressure_vector = pvt_data['pressure']
    
    # Determine the fluid key to use (default to 'brine' if the specified fluid is not found)
    fluid_key_to_use = fluid_key if fluid_key in pvt_data else 'brine'
    
    # Retrieve the density matrix for the specified fluid
    rho_matrix = pvt_data[fluid_key_to_use]['rho']
    
    # Create an interpolator for the density matrix
    rho_interpolator = RectBivariateSpline(pressure_vector, temperature_vector, rho_matrix)
    
    # Interpolate the density at the given pressure and temperature
    rho = rho_interpolator(pressure, temperature)[0, 0]
    
    return rho