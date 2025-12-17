import json
from pathlib import Path
from typing import Callable, Dict, Tuple, Union

import numpy as np
import pandas as pd
import scipy.constants as const
from scipy.integrate import solve_ivp
from scipy.interpolate import Akima1DInterpolator, RectBivariateSpline


def file_exists(files_path: list):
    """
    Receives a list of file paths and verify if the files in such
    path exists. Throws an exception in case the file is not found
    """
    for file in files_path:
        if not file.exists():
            raise FileNotFoundError(f"Required PVT file not found: {file}")


def get_mixture_info(pvt_path: Union[str, Path]):
    with open(f"{pvt_path}/metadata.json", "r") as file:
        mixture_info = json.load(file)
    return mixture_info


def load_envelopes(pvt_path: Union[str, Path]):
    bubble_point = np.loadtxt(pvt_path / "bubble_point.csv", delimiter=",", skiprows=1)
    dew_point = np.loadtxt(pvt_path / "dew_point.csv", delimiter=",", skiprows=1)

    return {"bubble_point": bubble_point, "dew_point": dew_point}


def set_envelope_interpolator(envelope_data: np.ndarray) -> Callable:
    """
    Sets an interpolator function for the envelope data
    """
    if envelope_data is None:
        return None

    T_envelope = envelope_data[:, 0]
    p_envelope = envelope_data[:, 1]

    # sort by T
    sorted_indices = np.argsort(T_envelope)

    T_envelope = T_envelope[sorted_indices]
    p_envelope = p_envelope[sorted_indices]

    envelope_interpolator = Akima1DInterpolator(T_envelope, p_envelope, method="makima", extrapolate=False)

    return envelope_interpolator


def corr_rhobrine_LaliberteCopper(salinity: float, temperature: np.ndarray, pressure: np.ndarray, rho_h2o: np.ndarray) -> np.ndarray:
    """
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
    """

    ## Laliberté and Cooper model for NaCl solutions
    # Laliberté and Cooper model: constants for NaCl
    c0 = -0.00433
    c1 = 0.06471
    c2 = 1.0166
    c3 = 0.014624
    c4 = 3315.6

    # NaCl concentration
    w = salinity / 100

    # Laliberté and Cooper model: Apparent density
    rho_app = (c0 * w + c1) * np.exp(0.000001 * (temperature + c4) ** 2) / (w + c2 + c3 * temperature)

    # Laliberté and Cooper model: Brine density
    rho_brine = 1 / (((1 - w) / rho_h2o) + (w / rho_app))

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
    rho_brine = np.loadtxt(brine_rho_file, delimiter=",")  # Assuming comma separation

    # Load metadata for brine
    brine_metadata = get_mixture_info(brine_path)

    # load brine envelopes if available
    brine_envelopes = load_envelopes(brine_path)

    bubble_point_brine = brine_envelopes.get("bubble_point", None)
    dew_point_brine = brine_envelopes.get("dew_point", None)

    bubble_point_brine_interpolator = set_envelope_interpolator(bubble_point_brine)
    dew_point_brine_interpolator = set_envelope_interpolator(dew_point_brine)

    # Create a dictionary for brine PVT data
    pvt_data = {
        "temperature": temperature_data,
        "pressure": pressure_data,
        "brine": {
            "rho": rho_brine,
            "metadata": brine_metadata,
            "bubble_point": bubble_point_brine_interpolator,
            "dew_point": dew_point_brine_interpolator,
        },
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
        rho_fluid = np.loadtxt(fluid_rho_file, delimiter=",")  # Assuming comma separation

        # Load metadata for the specified fluid
        fluid_metadata = get_mixture_info(fluid_path)

        # Load fluid envelopes if available
        fluid_envelopes = load_envelopes(fluid_path)
        fluid_bubble_point = fluid_envelopes.get("bubble_point", None)
        fluid_dew_point = fluid_envelopes.get("dew_point", None)

        fluid_bubble_point_interpolator = set_envelope_interpolator(fluid_bubble_point)
        fluid_dew_point_interpolator = set_envelope_interpolator(fluid_dew_point)

        # Add fluid data to the PVT data dictionary
        pvt_data[fluid_type] = {
            "rho": rho_fluid,
            "metadata": fluid_metadata,
            "bubble_point": fluid_bubble_point_interpolator,
            "dew_point": fluid_dew_point_interpolator,
        }

    return pvt_data


def _check_phase(P: float, T: float, T_crit: float, bubble_interp: Callable, dew_interp: Callable):
    """
    Return one of: "supercritical", "liquid", "gas", "two-phase", "unknown".
    Safely handles missing interpolators or NaN envelope values.
    """
    if (T_crit is not None) and (T >= T_crit):
        return "supercritical"

    # Evaluate envelopes
    p_bubble = np.nan
    p_dew = np.nan

    if bubble_interp is not None:
        try:
            p_bubble = float(bubble_interp(T))
        except Exception:
            p_bubble = np.nan

    if dew_interp is not None:
        try:
            p_dew = float(dew_interp(T))
        except Exception:
            p_dew = np.nan

    # If both envelopes invalid, return "unknown"
    if np.isnan(p_bubble) and np.isnan(p_dew):
        return "unknown"

    # Check liquid/gas/two-phase (prefer bubble/dew where available)
    if not np.isnan(p_bubble) and P > p_bubble:
        return "liquid"
    if not np.isnan(p_dew) and P < p_dew:
        return "gas"
    return "two-phase"


def _flag_phase_changes(P: float, phase: str, p_dew: float, p_bubble: float, tol: float = 0.665):
    """
    Return a flag string when approaching envelope within tol (bar),
    or None otherwise.
    """

    # Ensure p_dew / p_bubble may be NaN
    if phase == "gas" and not np.isnan(p_dew):
        diff_dew = p_dew - P
        if diff_dew <= tol:
            return "approaching dew point"
    if phase == "liquid" and not np.isnan(p_bubble):
        diff_bubble = P - p_bubble
        if diff_bubble <= tol:
            return "approaching bubble point"
    return None


def _eval_envelopes_for_temps(temps: np.ndarray, dew_interp: Callable, bub_interp: Callable, T_crit: float = None):
    """
    Returns (p_dew_vals, p_bub_vals) arrays aligned with temps.
    If an interpolator is None or fails or returns NaN, corresponding entries are NaN.
    Also mark values above T_crit as NaN.
    """
    # dew
    if dew_interp is not None:
        try:
            p_dew_vals = np.asarray(dew_interp(temps), dtype=float)
            # ensure shape matches temps
            if p_dew_vals.shape != temps.shape:
                p_dew_vals = np.broadcast_to(p_dew_vals, temps.shape).copy()
        except Exception:
            p_dew_vals = np.full(temps.shape, np.nan, dtype=float)
    else:
        p_dew_vals = np.full(temps.shape, np.nan, dtype=float)

    # bubble
    if bub_interp is not None:
        try:
            p_bub_vals = np.asarray(bub_interp(temps), dtype=float)
            if p_bub_vals.shape != temps.shape:
                p_bub_vals = np.broadcast_to(p_bub_vals, temps.shape).copy()
        except Exception:
            p_bub_vals = np.full(temps.shape, np.nan, dtype=float)
    else:
        p_bub_vals = np.full(temps.shape, np.nan, dtype=float)

    # invalidate above critical temperature if T_crit provided
    if T_crit is not None:
        invalid = temps >= T_crit
        p_dew_vals[invalid] = np.nan
        p_bub_vals[invalid] = np.nan

    return p_dew_vals, p_bub_vals


# Define _Pdz_odesys outside with necessary arguments
def _Pdz_odesys(z, y, depth_array: np.ndarray, temperature_array: np.ndarray, rho_interpolator: RectBivariateSpline, fluid_key: str) -> Tuple[float]:
    P = float(y[0])
    T = float(np.interp(z, depth_array, temperature_array))
    rho = get_rho_from_pvt_data(P, T, rho_interpolator)
    dPdz = rho * const.g / const.bar
    return np.array([dPdz])


def _integrate_pressure(
    init_curves: pd.DataFrame,
    reference_depth: float,
    reference_pressure: float,
    pvt_data: dict,
    fluid_key: str,
    interpolator: RectBivariateSpline,
    top_limit: float = None,
    bottom_limit: float = None,
) -> np.ndarray:
    """
    Simple integration to find pressure
    Starting point: reference_depth at reference pressure.  Reference temperature is found in init_curves
                    Then iterate upwards (up) or downwards (down) to subtract or add pressure.
                    Recalculates rho at each step.
    """
    init_curves = init_curves.copy()

    warnings = []
    warning_item = {"p": np.nan, "T": np.nan, "z": np.nan, "message": ""}

    depth_array = init_curves["depth"].values

    # metadata and envelopes
    T_crit = pvt_data[fluid_key]["metadata"].get("T_crit", None)
    bubble_point_interpolator = pvt_data[fluid_key].get("bubble_point", None)
    dew_point_interpolator = pvt_data[fluid_key].get("dew_point", None)

    # evaluate envelopes into DataFrame columns safely
    temps = init_curves["temperature"].values
    p_dew_vals, p_bub_vals = _eval_envelopes_for_temps(temps, dew_point_interpolator, bubble_point_interpolator, T_crit)
    init_curves["dew_point"] = p_dew_vals
    init_curves["bubble_point"] = p_bub_vals

    ## Initialization
    # New columns needed in DataFrame
    if fluid_key == "brine":
        fluid_label = fluid_key
    else:
        fluid_label = "fluid"

    colname_p = fluid_label + "_pressure"

    if colname_p not in init_curves.columns:
        init_curves[colname_p] = np.nan

    # Check and update limits
    if top_limit is None:
        top_limit = float(depth_array.min())

    if bottom_limit is None:
        bottom_limit = float(depth_array.max())

    # Check if reference_depth is within limits
    if reference_depth < top_limit or reference_depth > bottom_limit:
        raise ValueError("reference_depth is out of range")

    # Check if the reference_depth value is equal to bottom AND top limit
    if reference_depth == top_limit and reference_depth == bottom_limit:
        raise ValueError("reference_depth cannot be equal to bottom and top limit")

    # Check and integrate upwards if needed
    if reference_depth > top_limit:
        ivp_data = build_ivp_data(reference_depth, top_limit, reference_pressure, init_curves, pvt_data[fluid_key])
        ivp_solution = solve_ivp_with_data(ivp_data, interpolator, fluid_key)
        init_curves.loc[ivp_data["index"], colname_p] = ivp_solution.y[0]

        reference_temperature = np.interp(reference_depth, init_curves["depth"].values, init_curves["temperature"].values)

        reference_phase = _check_phase(reference_pressure, reference_temperature, T_crit, bubble_point_interpolator, dew_point_interpolator)

        init_curves["phase"] = init_curves.apply(
            lambda row: _check_phase(
                row[colname_p],
                row["temperature"],
                T_crit,
                bubble_point_interpolator,
                dew_point_interpolator,
            ),
            axis=1,
        )

        init_curves["phase_change_flag"] = init_curves.apply(
            lambda row: _flag_phase_changes(
                row[colname_p],
                row["phase"],
                row["dew_point"],
                row["bubble_point"],
            ),
            axis=1,
        )

        # Logging and trigger warnings
        supercritical_log = init_curves.loc[ivp_data["index"]][["phase"]].drop_duplicates()
        if reference_phase == "supercritical" and supercritical_log.shape[0] >= 2 and fluid_key != "brine":
            p_log = init_curves.loc[supercritical_log.index[1], colname_p].astype(float)
            T_log = init_curves.loc[supercritical_log.index[1], "temperature"].astype(float)
            z_log = init_curves.loc[supercritical_log.index[1], "depth"].astype(float)
            initial_phase = supercritical_log["phase"].values[0]
            final_phase = supercritical_log["phase"].values[1]
            message = f"Message: Fluid changed phase from {initial_phase} to {final_phase} at {z_log} mTVDMSL (P = {p_log:.1f} bar)."

            warnings.append({"p": p_log, "T": T_log, "z": z_log, "message": message})

        phase_log = init_curves.loc[ivp_data["index"]][["phase_change_flag"]].drop_duplicates().dropna()
        if phase_log.shape[0] >= 1 and fluid_key != "brine":
            p_log = init_curves.loc[phase_log.index[0], colname_p].astype(float)
            T_log = init_curves.loc[phase_log.index[0], "temperature"].astype(float)
            z_log = init_curves.loc[phase_log.index[0], "depth"].astype(float)
            phase = init_curves.loc[phase_log.index[0], "phase"]
            trigger = phase_log["phase_change_flag"].values[0]
            message = f"Warning: Fluid is in {phase} phase and is {trigger} at {z_log} m (P = {p_log:.1f} bar)."
            warnings.append({"p": p_log, "T": T_log, "z": z_log, "message": message})

    # Check and integrate downwards if needed
    if reference_depth <= bottom_limit:
        ivp_data = build_ivp_data(reference_depth, bottom_limit, reference_pressure, init_curves, pvt_data[fluid_key])
        ivp_solution = solve_ivp_with_data(ivp_data, interpolator, fluid_key)
        init_curves.loc[ivp_data["index"], colname_p] = ivp_solution.y[0]

    if reference_depth in init_curves["depth"].values:
        init_curves.loc[init_curves["depth"] == reference_depth, colname_p] = reference_pressure

    pressures = init_curves[colname_p].to_numpy(dtype=float)

    return pressures, warnings


def build_ivp_data(interval_start: float, interval_end: float, initial_state: float, init_curves: pd.DataFrame, pvt_data: dict) -> Dict:
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
    depth = init_curves["depth"].to_numpy()
    temperature = init_curves["temperature"].to_numpy()

    query = None
    if interval_start >= interval_end:
        query = init_curves[init_curves.depth < interval_start]
        query = query.sort_values("depth", ascending=False)
    elif interval_start <= interval_end:
        query = init_curves[init_curves.depth > interval_start]

    index = query.index.to_numpy()
    t_eval = query["depth"].to_numpy()

    ivp_data = {
        "interval_start": interval_start,
        "interval_end": interval_end,
        "initial_state": initial_state,
        "index": index,
        "t_eval": t_eval,
        "depth_array": depth,
        "temperature_array": temperature,
    }

    return ivp_data


def solve_ivp_with_data(ivp_data: dict, interpolator: RectBivariateSpline, fluid_key: str):
    ivp_solution = solve_ivp(
        fun=_Pdz_odesys,
        t_span=[ivp_data["interval_start"], ivp_data["interval_end"]],
        y0=[ivp_data["initial_state"]],
        t_eval=ivp_data["t_eval"],
        dense_output=True,
        args=(ivp_data["depth_array"], ivp_data["temperature_array"], interpolator, fluid_key),
        method="Radau",
    )
    return ivp_solution


def get_rho_from_pvt_data(pressure: float, temperature: float, rho_interpolator: RectBivariateSpline) -> float:
    # Interpolate the density at the given pressure and temperature
    rho = rho_interpolator(pressure, temperature)[0, 0]
    return rho
