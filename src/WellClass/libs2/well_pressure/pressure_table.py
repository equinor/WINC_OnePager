from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from scipy import constants as const
from scipy.interpolate import interp1d

from .pressure_utils import shmin_data_interpolator

if TYPE_CHECKING:
    from .fluid_pressure import FluidPressure


@dataclass
class PressureTable:
    name: str
    top_depth: float
    bottom_depth: float

    # Input parameters
    ground_elevation: float  # ground elevation (m)
    ground_temperature: float  # ground temperature (°C)
    geothermal_gradient: float  # geothermal gradient (°C/km)

    step: float = 1.0  # step size for depth array (m)
    rho_brine: float = 1030  # brine density (kg/m³)

    # Minimum horizontal stress parameters
    shmin_gradient: float = 0.1695
    shmin_data: list[list[float]] = field(default=None)  # depth values for shmin data

    # Computed arrays
    depth: np.ndarray = field(init=False)  # depth array (m)
    temperature: np.ndarray = field(init=False)  # temperature array (°C)
    hydrostatic_pressure: np.ndarray = field(init=False)  # hydrostatic pressure (bar)
    min_horizontal_stress: np.ndarray = field(init=False)  # minimum horizontal stress (bar)

    # Collection of fluid pressure scenarios
    scenarios: list[FluidPressure] = field(init=False, default_factory=list)

    def __post_init__(self) -> None:
        # Generate depth array
        self.depth = np.linspace(self.top_depth, self.bottom_depth, int((self.bottom_depth - self.top_depth) / self.step) + 1)

        self.temperature = self.compute_temperature()
        self.hydrostatic_pressure = self.compute_hydrostatic_pressure()
        self.min_horizontal_stress = self.compute_shmin()

    def compute_temperature(self) -> np.ndarray:
        """Compute temperature at a given depth using the geothermal gradient."""
        temp_array = self.ground_temperature + (self.geothermal_gradient * (self.depth - self.ground_elevation) / 1000)

        temp_array[self.depth < self.ground_elevation] = self.ground_temperature

        return temp_array

    def compute_hydrostatic_pressure(self) -> np.ndarray:
        pressure_pa = const.atm + self.depth * self.rho_brine * const.g
        return pressure_pa / 1e5

    def compute_shmin(self) -> np.ndarray:
        """Compute minimum horizontal stress at a given depth."""
        ground_pressure = np.interp(self.ground_elevation, self.depth, self.hydrostatic_pressure)

        if self.shmin_data is not None:
            return shmin_data_interpolator(
                shmin_data=self.shmin_data,
                depth_array=self.depth,
                ground_elevation=self.ground_elevation,
                ground_pressure_bar=ground_pressure,
                surface_pressure_bar=const.atm / 1e5,
            )

        if self.shmin_gradient is not None:
            shmin = ground_pressure + self.shmin_gradient * (self.depth - self.ground_elevation)
            shmin[self.depth < self.ground_elevation] = self.hydrostatic_pressure[self.depth < self.ground_elevation]
            return shmin

        raise ValueError("Either shmin_gradient or shmin_data must be provided.")

    def get_values_at_depth(self, depth_value: float) -> dict:
        """Interpolate all curves at a given depth."""

        def interp(arr: np.ndarray) -> float:
            f = interp1d(self.depth, arr, bounds_error=False, fill_value="extrapolate")
            return float(f(depth_value))

        return {
            "temperature": interp(self.temperature),
            "hydrostatic_pressure": interp(self.hydrostatic_pressure),
            "min_horizontal_stress": interp(self.min_horizontal_stress),
        }

    def add_scenario(self, **kwargs: float | str | None) -> FluidPressure:
        """Create a FluidPressure scenario, register it, and return it."""
        from .fluid_pressure import FluidPressure  # noqa: PLC0415

        scenario = FluidPressure(table=self, **kwargs)
        self.scenarios.append(scenario)
        return scenario
