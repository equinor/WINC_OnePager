"""
WellConditions: Manages the initial physical conditions of a well.

This class is responsible for computing depth, temperature, hydrostatic pressure,
and minimum horizontal stress curves for a well based on its geometry and thermal properties.
"""

import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd
import scipy.constants as const
from scipy.interpolate import interp1d

from ..pvt.pvt import _integrate_pressure

# Constants
SHMIN_FAC = 0.1695  # empirical factor for Shmin calculation


@dataclass
class WellConditions:
    """
    Computes and manages initial physical condition curves for a well.

    Attributes:
        sf_depth_msl (float): Seafloor depth relative to mean sea level (MSL) in meters.
        well_td_rkb (float): Total depth of the well relative to the rotary kelly bushing (RKB) in meters.
        well_rkb (float): Rotary kelly bushing elevation in meters.
        sf_temp (float): Seafloor temperature in degrees Celsius.
        geo_tgrad (float): Geothermal gradient in degrees Celsius per kilometer.
        pvt_data (dict): PVT data dictionary containing pressure, temperature, and brine data.
        brine_interpolator: Interpolator for brine density.
        rho_brine (float, optional): Constant brine density if provided.
        salinity (float): Salinity percentage (default 3.5%).
        shmin_gradient (float): Gradient for Shmin calculation.
        ip_shmin_data (np.ndarray, optional): User-provided Shmin data.
        curves (pd.DataFrame): Computed initial curves.
    """

    sf_depth_msl: float
    well_td_rkb: float
    well_rkb: float
    sf_temp: float
    geo_tgrad: float
    pvt_data: dict
    brine_interpolator: object
    rho_brine: Optional[float] = None
    salinity: float = 3.5
    shmin_gradient: float = SHMIN_FAC
    ip_shmin_data: Optional[np.ndarray] = None
    curves: pd.DataFrame = field(init=False)

    def __post_init__(self):
        """Compute all initial condition curves after instantiation."""
        self._curves = self._compute_all_curves()
        # Add validation
        assert not self._curves['depth'].isna().any(), "Depth curve has NaN values"
        assert not self._curves['temperature'].isna().any(), "Temperature curve has NaN values"

    def _compute_all_curves(self) -> pd.DataFrame:
        """
        Computes all initial curves: depth, temperature, hydrostatic pressure, and Shmin.

        Returns:
            pd.DataFrame: DataFrame containing all computed curves.
        """
        depth_curve = self._calculate_depth_curve()
        curves = pd.DataFrame(
            {
                "depth": depth_curve,
                "temperature": np.nan,
                "hydrostatic_pressure": np.nan,
                "min_horizontal_stress": np.nan,
            }
        )

        curves["temperature"] = self._calculate_temperature_curve(depth_curve)
        curves["hydrostatic_pressure"] = self._calculate_hydrostatic_pressure(curves)
        curves["min_horizontal_stress"] = self._calculate_shmin(depth_curve, curves["hydrostatic_pressure"])

        return curves

    def _calculate_depth_curve(self) -> np.ndarray:
        """
        Calculate depth vector from MSL downwards.

        Returns:
            np.ndarray: Array of depth values.
        """
        dz = 1.0
        td_msl = self.well_td_rkb - self.well_rkb
        z_bottom = int(td_msl) + 500
        z_vec = np.arange(0, z_bottom, dz)
        return z_vec

    def _calculate_temperature_curve(self, depth_curve: np.ndarray) -> np.ndarray:
        """
        Calculate temperature curve based on geothermal gradient.

        Args:
            depth_curve: Array of depth values.

        Returns:
            np.ndarray: Array of temperature values.
        """
        temperature_curve = self.sf_temp + np.maximum(0, (self.geo_tgrad * (depth_curve - self.sf_depth_msl)) / 1e3)
        return temperature_curve

    def _calculate_hydrostatic_pressure(self, curves: pd.DataFrame) -> np.ndarray:
        """
        Calculate hydrostatic pressure curve.

        Args:
            curves: DataFrame containing depth and temperature curves.

        Returns:
            np.ndarray: Array of hydrostatic pressure values.
        """
        if self.rho_brine is not None:
            # Use constant density
            hydrostatic_pressure = (const.atm + curves["depth"] * self.rho_brine * const.g) / const.bar
        else:
            # Use PVT data integration
            hydrostatic_pressure, _ = _integrate_pressure(
                init_curves=curves,
                reference_depth=0,
                reference_pressure=const.atm / const.bar,
                pvt_data=self.pvt_data,
                fluid_key="brine",
                interpolator=self.brine_interpolator,
            )

        return hydrostatic_pressure

    def _calculate_shmin(self, depth_array: np.ndarray, hydrostatic_pressure_curve: np.ndarray) -> np.ndarray:
        """
        Compute minimum horizontal stress (Shmin) based on user data or empirical formula.

        Args:
            depth_array: Array of depth values.
            hydrostatic_pressure_curve: Array of hydrostatic pressure values.

        Returns:
            np.ndarray: Array of Shmin values.
        """
        pressure_ml = np.interp(self.sf_depth_msl, depth_array, hydrostatic_pressure_curve)

        if self.ip_shmin_data is not None:
            shmin_curve = self._interpolate_user_shmin(depth_array, hydrostatic_pressure_curve, pressure_ml)
        else:
            shmin_curve = self._calculate_empirical_shmin(depth_array, pressure_ml)

        # Ensure Shmin is not below hydrostatic pressure at any depth
        shmin_curve[depth_array < self.sf_depth_msl] = hydrostatic_pressure_curve[depth_array < self.sf_depth_msl]

        return shmin_curve

    def _interpolate_user_shmin(self, depth_array: np.ndarray, hydrostatic_pressure_curve: np.ndarray, pressure_ml: float) -> np.ndarray:
        """
        Interpolate user-provided Shmin data.

        Args:
            depth_array: Array of depth values.
            hydrostatic_pressure_curve: Array of hydrostatic pressure values.
            pressure_ml: Pressure at mudline.

        Returns:
            np.ndarray: Interpolated Shmin values.
        """
        depth_values, shmin_values = self.ip_shmin_data.T

        if min(depth_values) > self.sf_depth_msl:
            warnings.warn(
                f"No Shmin data between seafloor depth ({self.sf_depth_msl}) and minimum provided depth ({min(depth_values)}). "
                "Extrapolating using hydrostatic pressure at seafloor."
            )
        else:
            # Filter to depths >= seafloor
            mask = depth_values >= self.sf_depth_msl
            depth_values = depth_values[mask]
            shmin_values = shmin_values[mask]

        # Remove duplicates
        unique_indices = np.unique(depth_values, return_index=True)[1]
        depth_values = depth_values[np.sort(unique_indices)]
        shmin_values = shmin_values[np.sort(unique_indices)]

        # Insert seafloor and MSL values
        depth_values = np.insert(depth_values, 0, [0, self.sf_depth_msl])
        shmin_values = np.insert(shmin_values, 0, [hydrostatic_pressure_curve[0], pressure_ml])

        shmin_interpolator = interp1d(depth_values, shmin_values, bounds_error=False, fill_value="extrapolate")
        shmin_curve = shmin_interpolator(depth_array)

        # Check if Shmin values are problematically below hydrostatic
        below_sf = depth_array > self.sf_depth_msl
        if below_sf.any():
            check_ratio = (shmin_curve[below_sf] < hydrostatic_pressure_curve[below_sf]).sum() / below_sf.sum()
            if np.isclose(check_ratio, 1.0):
                warnings.warn(
                    f"Shmin values below seafloor depth ({self.sf_depth_msl}) are below hydrostatic pressure. "
                    "Issues with PREDICT Shmin data. Shmin values will be set to hydrostatic pressure."
                )

        # Ensure Shmin >= hydrostatic everywhere
        shmin_curve[shmin_curve < hydrostatic_pressure_curve] = hydrostatic_pressure_curve[shmin_curve < hydrostatic_pressure_curve]

        return shmin_curve

    def _calculate_empirical_shmin(self, depth_array: np.ndarray, pressure_ml: float) -> np.ndarray:
        """
        Calculate Shmin using empirical formula.

        Args:
            depth_array: Array of depth values.
            pressure_ml: Pressure at mudline.

        Returns:
            np.ndarray: Empirically calculated Shmin values.
        """
        depth_ml = depth_array - self.sf_depth_msl
        shmin_curve = pressure_ml + depth_ml * self.shmin_gradient
        return shmin_curve

    def get_curves(self) -> pd.DataFrame:
        """
        Get a copy of the computed curves.

        Returns:
            pd.DataFrame: Copy of curves DataFrame.
        """
        return self.curves.copy(deep=True)
