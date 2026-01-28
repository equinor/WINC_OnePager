"""
PVTDataManager: Manages PVT (Pressure-Volume-Temperature) data loading and interpolation.

This class handles loading fluid property data and creating interpolators for
density calculations across pressure and temperature ranges.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, Union

import numpy as np
from scipy.interpolate import RectBivariateSpline

from ..pvt.pvt import corr_rhobrine_LaliberteCopper, load_pvt_data


@dataclass
class PVTDataManager:
    """
    Manages PVT data loading and interpolation for fluids.

    Attributes:
        pvt_path (Union[str, Path]): Path to PVT data files.
        fluid_type (str, optional): Type of fluid to load.
        specific_gravity (float, optional): Specific gravity if using constant density.
        salinity (float): Salinity percentage for brine correction (default 3.5%).
        pvt_data (dict): Loaded PVT data.
        brine_interpolator: Interpolator for brine density.
        fluid_interpolator: Interpolator for fluid density (None if using specific_gravity).
        fluid_composition (str, optional): Composition of the fluid.
    """

    pvt_path: Union[str, Path]
    fluid_type: Optional[str] = None
    specific_gravity: Optional[float] = None
    salinity: float = 3.5
    pvt_data: Dict[str, Dict[str, np.ndarray]] = field(init=False)
    brine_interpolator: Optional[RectBivariateSpline] = field(init=False)
    fluid_interpolator: Optional[RectBivariateSpline] = field(init=False)
    fluid_composition: Optional[str] = field(init=False, default=None)

    def __post_init__(self):
        """Load PVT data and create interpolators after instantiation."""
        self._validate_inputs()
        self.pvt_data = self._load_pvt_data()
        self._initialize_interpolators()

    def _validate_inputs(self):
        """Validate that either fluid_type or specific_gravity is provided, but not both."""
        if self.fluid_type is None and self.specific_gravity is None:
            raise ValueError("Either 'fluid_type' or 'specific_gravity' must be provided.")
        if self.fluid_type is not None and self.specific_gravity is not None:
            raise ValueError("Cannot specify both 'fluid_type' and 'specific_gravity'.")

    def _load_pvt_data(self) -> dict:
        """
        Load PVT data based on configuration.

        Returns:
            dict: PVT data dictionary.
        """
        if self.specific_gravity is not None:
            # Only load brine data when using specific gravity
            return load_pvt_data(self.pvt_path, load_fluid=False)
        else:
            # Load fluid and brine data
            pvt_data = load_pvt_data(self.pvt_path, self.fluid_type, load_fluid=True)
            self.fluid_composition = pvt_data[self.fluid_type]["metadata"]["composition"]
            return pvt_data

    def _initialize_interpolators(self):
        """Create interpolators for brine and fluid density."""
        pressure_vector = self.pvt_data["pressure"]
        temperature_vector = self.pvt_data["temperature"]

        # Create brine interpolator
        self._initialize_brine_interpolator(pressure_vector, temperature_vector)

        # Create fluid interpolator (if not using specific gravity)
        if self.specific_gravity is None:
            self._initialize_fluid_interpolator(pressure_vector, temperature_vector)
        else:
            self.fluid_interpolator = None

    def _initialize_brine_interpolator(self, pressure_vector: np.ndarray, temperature_vector: np.ndarray):
        """
        Initialize the brine density interpolator with salinity correction.

        Args:
            pressure_vector: Array of pressure values.
            temperature_vector: Array of temperature values.
        """
        water_rho_matrix = self.pvt_data["brine"]["rho"]

        # Create meshgrid for salinity correction
        t_matrix, p_matrix = np.meshgrid(pressure_vector, temperature_vector)

        # Apply Laliberté and Cooper salinity correction
        brine_rho_matrix = corr_rhobrine_LaliberteCopper(self.salinity, t_matrix, p_matrix, water_rho_matrix)

        self.brine_interpolator = RectBivariateSpline(pressure_vector, temperature_vector, brine_rho_matrix)

    def _initialize_fluid_interpolator(self, pressure_vector: np.ndarray, temperature_vector: np.ndarray):
        """
        Initialize the fluid density interpolator.

        Args:
            pressure_vector: Array of pressure values.
            temperature_vector: Array of temperature values.
        """
        fluid_rho_matrix = self.pvt_data[self.fluid_type]["rho"]
        self.fluid_interpolator = RectBivariateSpline(pressure_vector, temperature_vector, fluid_rho_matrix)

    def load_additional_fluid(self, fluid_type: str) -> tuple:
        """
        Load an additional fluid type and create its interpolator.

        Args:
            fluid_type: Name of the fluid type to load.

        Returns:
            tuple: (fluid_interpolator, fluid_composition, pvt_data_for_fluid)
        """
        pvt_data = load_pvt_data(self.pvt_path, fluid_type, load_fluid=True)

        pressure_vector = pvt_data["pressure"]
        temperature_vector = pvt_data["temperature"]
        rho_matrix = pvt_data[fluid_type]["rho"]

        fluid_interpolator = RectBivariateSpline(pressure_vector, temperature_vector, rho_matrix)
        fluid_composition = pvt_data[fluid_type]["metadata"]["composition"]

        return fluid_interpolator, fluid_composition, pvt_data

    def get_brine_interpolator(self) -> RectBivariateSpline:
        """
        Get the brine interpolator.

        Returns:
            RectBivariateSpline: Brine density interpolator.
        """
        return self.brine_interpolator

    def get_fluid_interpolator(self) -> Optional[RectBivariateSpline]:
        """
        Get the fluid interpolator.

        Returns:
            Optional[RectBivariateSpline]: Fluid density interpolator, or None if using specific_gravity.
        """
        return self.fluid_interpolator

    def get_pvt_data(self) -> dict:
        """
        Get the PVT data dictionary.

        Returns:
            dict: PVT data.
        """
        return self.pvt_data
