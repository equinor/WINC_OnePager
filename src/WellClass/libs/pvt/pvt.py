from typing import Union, Dict, Tuple, Callable
from pathlib import Path
import numpy as np
import pandas as pd
import json
from scipy.interpolate import RectBivariateSpline
from scipy.integrate import solve_ivp

import scipy.constants as const


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



    # Check that fluid-specific files exist before proceeding
    for file in [brine_metadata_file, brine_rho_file]:
        if not file.exists():
            raise FileNotFoundError(f"Required PVT file not found: {file}")

    # Load density data
    rho_brine = np.loadtxt(brine_rho_file, delimiter=',')  # Assuming comma separation

    # Load metadata for brine
    brine_metadata = get_mixture_info(brine_path)


    # Create a dictionary for brine PVT data
    pvt_data = {
        "temperature": t,
        "pressure": p,
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
        for file in [fluid_metadata_file, fluid_rho_file]:
            if not file.exists():
                raise FileNotFoundError(f"Required PVT file not found: {file}")

        # Load density data for the specified fluid
        rho_fluid = np.loadtxt(fluid_rho_file, delimiter=',')  # Assuming comma separation

        # Load metadata for the specified fluid
        fluid_metadata = get_mixture_info(fluid_path)

        # Add fluid data to the PVT data dictionary
        pvt_data[fluid_type] = { "rho": rho_fluid,
                                "metadata": fluid_metadata }


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
    return dPdz,

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
    
    
    # ip_arrays = {'depth': depth_array, 'temperature': temperature_array}
    # init_curves = pd.DataFrame(data = ip_arrays)

    depth_array = init_curves['depth'].values
    temperature_array = init_curves['temperature'].values


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
                         args=(depth_array, temperature_array, interpolator, fluid_key),
                         method = 'Radau')
    
        init_curves.loc[query_up.index, colname_p] = solution_up.y[0]

    # Check and integrate downwards if needed
    if reference_depth <= bottom_limit:

        query_down = init_curves[init_curves.depth > reference_depth]
        solution_down = solve_ivp(_Pdz_odesys, [reference_depth, bottom_limit], [reference_pressure], 
                         t_eval=query_down['depth'].values, dense_output=True,
                         args=(depth_array, temperature_array, interpolator, fluid_key),
                         method = 'Radau')

        init_curves.loc[query_down.index, colname_p] = solution_down.y[0]

    if reference_depth in init_curves['depth'].values:
        init_curves.loc[init_curves['depth'] == reference_depth, colname_p] = reference_pressure

    return init_curves[colname_p].values.astype(float)




def get_rho_from_pvt_data(pressure: float, temperature: float, rho_interpolator: RectBivariateSpline) -> float:
    # Interpolate the density at the given pressure and temperature
    rho = rho_interpolator(pressure, temperature)[0, 0]
    return rho