from dataclasses import dataclass, field
from typing import Union, Dict
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from .new_PressureScenarioManager import PressureScenarioManager
from .new_aux_func import load_pvt_data, compute_hydrostatic_pressure

# Constants
SHMIN_NAME = 'Shmin'
SHMIN_FAC = 0.1695  # empirical factor for Shmin calculation

@dataclass
class Pressure:
    sf_depth_msl: float
    well_td_rkb: float
    well_rkb: float
    sf_temp: float
    geo_tgrad: float  # Geothermal gradient in degC/km
    pvt_path: Union[str, Path]
    fluid_type: str
    z_fluid_contact: float = None
    p_fluid_contact: float = None
    z_resrv: float = None
    p_resrv: float = None
    pvt_data: Dict[str, Dict[str, np.ndarray]] = field(init=False)
    ip_shmin_data: np.ndarray = field(default=None) # User-provided Shmin data
    init_curves: pd.DataFrame = field(init=False)
    scenario_manager: PressureScenarioManager = field(init=False)
    input_scenarios: Dict[str, Union[float, str, None]] = field(default_factory=dict)
    default_hs_scenario: bool = True
    collated_profiles: pd.DataFrame = field(init=None)


    def __post_init__(self):
        self.pvt_data = load_pvt_data(self.pvt_path, self.fluid_type)
        self.init_curves = self._compute_init_curves()
        self.scenario_manager = PressureScenarioManager()
        self.manage_scenarios()

    def _get_pvt_data(self) -> Dict[str, np.ndarray]:
        # Retrieve PVT data only once and store it for use throughout the class
        # This method will call a function that loads the PVT data files and returns a dictionary
        return load_pvt_data(self.pvt_path)

    def _compute_init_curves(self) -> pd.DataFrame:
        # Logic to compute depth, temperature, hydrostatic pressure, and Shmin
        # Return a pandas DataFrame with these initial curves
        depth_curve = self._calculate_depth_curve()
        temperature_curve = self._calculate_temperature_curve(depth_curve)
        hydrostatic_pressure_curve = self._calculate_hydrostatic_pressure(depth_curve, temperature_curve)
        min_horizontal_stress = self._calculate_shmin(depth_curve, hydrostatic_pressure_curve)
        
        # Combine into a single DataFrame
        init_curves = pd.DataFrame({
            'depth': depth_curve,
            'temperature': temperature_curve,
            'hydrostatic_pressure': hydrostatic_pressure_curve,
            'min_horizontal_stress': min_horizontal_stress
        })

        return init_curves


    def _calculate_depth_curve(self) -> np.ndarray:
        #Make the depth-vector from msl and downwards
        dz = 1.0
        td_msl = self.well_td_rkb - self.well_rkb
        z_bottom = int(td_msl)+500
        z_vec  = np.arange(0, z_bottom, dz)

        return z_vec
    
    def _calculate_temperature_curve(self, depth_curve: np.ndarray) -> np.ndarray:
        # Calculate depth sample points based on well depth and seabed depth
        temperature_curve = self.sf_temp + np.maximum(0, (self.geo_tgrad * (depth_curve - self.sf_depth_msl)) / 1e3)
        
        return temperature_curve

    def _calculate_hydrostatic_pressure(self, depth_curve: np.ndarray, temperature_curve: np.ndarray) -> np.ndarray:
         # Compute hydrostatic pressure using density data from pvt_data and integrate pressure
         return compute_hydrostatic_pressure(depth_curve, temperature_curve, self.pvt_data)

    def _calculate_shmin(self, depth_array: np.ndarray, hydrostatic_pressure_curve: np.ndarray) -> np.ndarray:
        """
        Compute Shmin based on user-provided data or an empirical formula.

        Args:
            depth_array (np.ndarray): Array of depth values.
            hydrostatic_pressure_curve (np.ndarray): Array of hydrostatic pressure values.

        Returns:
            np.ndarray: Array of Shmin values.
        """
        if self.ip_shmin_data is not None:
            # User has provided custom Shmin data, so interpolate it
            depth_values, shmin_values = self.user_shmin_data.T
            shmin_interpolator = interp1d(depth_values, shmin_values, bounds_error=False, fill_value="extrapolate")
            shmin_curve = shmin_interpolator(depth_array)
        else:
            # No user data provided, use the empirical formula
            shmin_curve = np.copy(hydrostatic_pressure_curve) # Copy hydrostatic pressure as a starting point
            below_sf_mask = depth_array >= self.sf_depth_msl # Mask for depths below seafloor
            shmin_curve[below_sf_mask] += (depth_array[below_sf_mask] - self.sf_depth_msl) * SHMIN_FAC

        return shmin_curve

    def add_scenario(self, scenario_name: str, **kwargs):
        # Ensure init_curves is included in kwargs or add it if not present
        kwargs.setdefault('init_curves', self.init_curves)
        kwargs.setdefault('fluid_type', self.fluid_type)
        kwargs.setdefault('pvt_data', self.pvt_data)
        kwargs.setdefault('z_fluid_contact', self.z_fluid_contact)
        
        # Create a new PressureScenario instance with the given name and parameters
        scenario = self.scenario_manager.create_scenario(name=scenario_name, **kwargs)
        scenario.compute_pressure_profile()

    def manage_scenarios(self):
        # Check if the first scenario is 'None' and should be computed as hydrostatic
        if self.z_fluid_contact is not None and self.input_scenarios:
            # Skip 'depth_msl' and find the first actual scenario name and pressure
            scenarios_iter = ((name, pressure) for name, pressure in self.input_scenarios.items() if name.startswith('RP'))
            first_scenario_name, first_scenario_pressure = next(scenarios_iter, (None, None))
            
            # Check if the first scenario pressure is None and should be computed as hydrostatic
            if first_scenario_pressure is None:
                # Compute hydrostatic scenario named as the first scenario
                self.add_scenario(
                    scenario_name=first_scenario_name,
                    from_resrvr=True,
                    z_fluid_contact=self.z_fluid_contact,
                    p_delta=0  # Assuming hydrostatic pressure
                )

            
            # Parse the remaining scenarios and create PressureScenario instances
            for sc_name, sc_pressure in scenarios_iter:
                if sc_pressure is not None:
                    p_delta = self._parse_pressure_magnitude(sc_pressure)
                    self.add_scenario(
                        scenario_name=sc_name,
                        from_resrvr=True,
                        z_fluid_contact=self.z_fluid_contact,
                        p_delta=p_delta
                    )

        
        elif self.z_fluid_contact is not None:
            # Handle default hydrostatic scenario if no input scenarios are provided
            self.add_scenario(
                scenario_name='hydrostatic',
                from_resrvr=True,
                z_fluid_contact=self.z_fluid_contact
            )
            
    def collate_all_profiles(self):
        # Call the collate_scenario_profiles method from the scenario manager
        common_data = self.init_curves[['depth', 'temperature', 'hydrostatic_pressure', 'min_horizontal_stress']]
        self.collated_profiles = self.scenario_manager.collate_scenario_profiles(common_data)

    
    def _parse_pressure_magnitude(self, sc_pressure: Union[str, float]) -> float:
        if isinstance(sc_pressure, str):
            # Convert string to float, handling any special string cases here
            sc_pressure = float(sc_pressure)

        return sc_pressure