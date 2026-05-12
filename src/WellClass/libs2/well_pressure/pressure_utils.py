import warnings

import numpy as np
from scipy.interpolate import interp1d


def shmin_data_interpolator(
    shmin_data: list[list[float]], depth_array: np.ndarray, ground_elevation: float, ground_pressure_bar: float, surface_pressure_bar: float
) -> np.ndarray:
    shmin_data = np.array(shmin_data)
    depth_values, shmin_values = shmin_data.T

    # Remove duplicates while keeping the first occurrence
    unique_indices = np.unique(depth_values, return_index=True)[1]
    depth_values = depth_values[np.sort(unique_indices)]
    shmin_values = shmin_values[np.sort(unique_indices)]

    if min(depth_values) > ground_elevation:
        warnings.warn(
            f"No Shmin data between seafloor depth ({ground_elevation}) and minimum provided depth ({min(depth_values)}). "
            "Extrapolating using hydrostatic pressure at seafloor."
        )

    else:
        filtered_depth_values = depth_values[depth_values >= ground_elevation]
        filtered_shmin_values = shmin_values[depth_values >= ground_elevation]

        depth_values = filtered_depth_values
        shmin_values = filtered_shmin_values

    # Insert the seafloor depth and pressure at mudline into the arrays
    depth_values = np.insert(depth_values, 0, ground_elevation)
    shmin_values = np.insert(shmin_values, 0, ground_pressure_bar)

    # Insert values af MSL
    depth_values = np.insert(depth_values, 0, 0)
    shmin_values = np.insert(shmin_values, 0, surface_pressure_bar)

    shmin_interpolator = interp1d(depth_values, shmin_values, bounds_error=False, fill_value="extrapolate")

    return shmin_interpolator(depth_array)
