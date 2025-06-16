from dataclasses import dataclass, field
from typing import Union, Dict, Optional
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d, RectBivariateSpline

from ..well_class.well_class import Well


from .PressureScenarioManager import PressureScenarioManager
from ..pvt.pvt import (load_pvt_data, 
                          compute_hydrostatic_pressure, 
                          get_rho_from_pvt_data,
                          corr_rhobrine_LaliberteCopper)
from .barrier_pressure import leakage_proxy
import warnings

# Constants
SHMIN_NAME = 'Shmin'
SHMIN_FAC = 0.1695  # empirical factor for Shmin calculation

@dataclass
class Pressure:
    """
    Manages pressure scenarios for a legacy well based on geothermal gradient and fluid properties.

    Attributes:
        sf_depth_msl (float): Seafloor depth relative to mean sea level (MSL) in meters.
        well_td_rkb (float): Total depth of the well relative to the rotary kelly bushing (RKB) in meters.
        well_rkb (float): Rotary kelly bushing elevation in meters.
        sf_temp (float): Seafloor temperature in degrees Celsius.
        geo_tgrad (float): Geothermal gradient in degrees Celsius per kilometer.
        pvt_path (Union[str, Path]): Path to the PVT data files.
        fluid_type (str): Name of the fluid stored in PVT data to use in pressure scenarios.
        z_fluid_contact (float, optional): Depth of the fluid contact in meters.
        p_fluid_contact (float, optional): Pressure at the fluid contact in bars.
        z_resrv (float, optional): Depth of the reservoir in meters.
        p_resrv (float, optional): Pressure at the reservoir in bars.
        fluid_composition (str, optional): Composition of the fluid.
        pvt_data (Dict[str, Dict[str, np.ndarray]]): PVT data for the fluid and brine.
        brine_interpolator (RectBivariateSpline): Interpolator for brine PVT data.
        fluid_interpolator (RectBivariateSpline): Interpolator for fluid PVT data.
        ip_shmin_data (np.ndarray, optional): User-provided minimum horizontal stress data.
        init_curves (pd.DataFrame): Initial curves for depth, temperature, hydrostatic pressure, and minimum horizontal stress.
        scenario_manager (PressureScenarioManager): Manager for pressure scenarios.
        input_scenarios (Dict[str, Union[float, str, None]]): Input scenarios for pressure calculations.
        default_hs_scenario (bool): Flag indicating whether to compute the default hydrostatic scenario.
        salinity (float): Salinity of the fluid in percentage (default is 3.5% for seawater).
        shmin_gradient (float): Gradient for Shmin calculation (default is SHMIN_FAC).
    """


    sf_depth_msl: float
    well_td_rkb: float
    well_rkb: float
    sf_temp: float
    geo_tgrad: float  # Geothermal gradient in degC/km
    pvt_path: Union[str, Path]
    fluid_type: str
    z_fluid_contact: Optional[float] = None
    p_fluid_contact: Optional[float] = None
    z_resrv: Optional[float] = None
    p_resrv: Optional[float] = None
    fluid_composition: Optional[str] = None
    specific_gravity: Optional[float] = None
    pvt_data: Dict[str, Dict[str, np.ndarray]] = field(init=False)
    brine_interpolator: RectBivariateSpline = field(init=False)
    fluid_interpolator: RectBivariateSpline = field(init=False)
    ip_shmin_data: np.ndarray = field(default=None) # User-provided Shmin data
    init_curves: pd.DataFrame = field(init=False)
    scenario_manager: PressureScenarioManager = field(init=False)
    input_scenarios: Dict[str, Optional[Union[float, str]]] = field(default_factory=dict)
    default_hs_scenario: bool = True
    salinity: float = 3.5  # Salinity in percentage (default is 3.5% for seawater)
    shmin_gradient: float = SHMIN_FAC  # Gradient for Shmin calculation


    def __post_init__(self):
        """
        Initializes PVT data, interpolators, initial curves, and manages scenarios 
        after class instantiation.
        """
        self.pvt_data = self._load_pvt_data()
        self._initialize_interpolators()
        self.init_curves = self._compute_init_curves()
        self.scenario_manager = PressureScenarioManager()
        self.manage_scenarios()

    def _load_pvt_data(self):
        """Load PVT data based on specific gravity and fluid type."""
        if self.specific_gravity is not None:
            return load_pvt_data(self.pvt_path, load_fluid=False)
        else:
            pvt_data = load_pvt_data(self.pvt_path, self.fluid_type, load_fluid=True)
            self.fluid_composition = pvt_data[self.fluid_type]['metadata']['composition']
            return pvt_data

    def _compute_init_curves(self) -> pd.DataFrame:
        """
        Computes initial curves for depth, temperature, hydrostatic pressure, and minimum horizontal stress.

        Returns:
            pd.DataFrame: DataFrame containing the initial curves.
        """
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

    def _initialize_interpolators(self, load_fluid: bool = True):
        # Initialize the brine interpolator
        pressure_vector = self.pvt_data['pressure']
        temperature_vector = self.pvt_data['temperature']

        water_rho_matrix = self.pvt_data['brine']['rho']

        t_matrix, p_matrix = np.meshgrid(pressure_vector, temperature_vector)

        # Correct water density for salinity using Laliberté and Cooper model
        brine_rho_matrix = corr_rhobrine_LaliberteCopper(self.salinity, t_matrix, p_matrix, water_rho_matrix)

        self.brine_interpolator = RectBivariateSpline(pressure_vector, temperature_vector, brine_rho_matrix)

        # Initialize the fluid interpolator based on the fluid type
        if self.specific_gravity is None:
            fluid_rho_matrix = self.pvt_data[self.fluid_type]['rho']
            self.fluid_interpolator = RectBivariateSpline(pressure_vector, temperature_vector, fluid_rho_matrix)
        else:
            # If specific gravity is provided, set fluid_interpolator to None
            self.fluid_interpolator = None


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
        pressure_ml = np.interp(self.sf_depth_msl, depth_array,  hydrostatic_pressure_curve)


        if self.ip_shmin_data is not None:
            # User has provided custom Shmin data, so interpolate it
            depth_values, shmin_values = self.ip_shmin_data.T

            if min(depth_values) > self.sf_depth_msl:
                print(min(depth_values))
                warnings.warn(
                    f"No Shmin data between seafloor depth ({self.sf_depth_msl}) and minimum provided depth ({min(depth_values)}). "
                    "Extrapolating using hydrostatic pressure at seafloor."
                )
                depth_values = np.insert(depth_values, 0, self.sf_depth_msl)
                shmin_values = np.insert(shmin_values, 0, pressure_ml)

            

            else:
                depth_values = depth_values[depth_values >= self.sf_depth_msl]
                shmin_values = shmin_values[depth_values >= self.sf_depth_msl]
                

            shmin_interpolator = interp1d(depth_values, shmin_values, bounds_error=False, fill_value="extrapolate")


            shmin_curve = shmin_interpolator(depth_array)

            check_shmin_values = (shmin_curve[depth_array>self.sf_depth_msl] < hydrostatic_pressure_curve[depth_array>self.sf_depth_msl]).sum() / len(shmin_curve[depth_array>self.sf_depth_msl])

            if np.isclose(check_shmin_values, 1.0):
                warnings.warn(
                    f"Shmin values below seafloor depth ({self.sf_depth_msl}) are below hydrostatic pressure.\nIssues with PREDICT Shmin data. Shmin values will be set to hydrostatic pressure."
                )

            # Ensure Shmin is equual to hydrostatic pressure above seafloor depth
            shmin_curve[shmin_curve<hydrostatic_pressure_curve] = hydrostatic_pressure_curve[shmin_curve<hydrostatic_pressure_curve]  
            
            



        else:
            # No user data provided, use the empirical formula
            
            # interpolate hydrostatic pressure at mudline depth (seafloor depth)
            pressure_ml = np.interp(self.sf_depth_msl, depth_array,  hydrostatic_pressure_curve)
            
            depth_ml = depth_array - self.sf_depth_msl  # depth below mean sea level

            shmin_curve = pressure_ml + depth_ml*self.shmin_gradient
            # shmin_curve[shmin_curve<0] = hydrostatic_pressure_curve[shmin_curve<0]

        # Ensure Shmin is not below hydrostatic pressure at any depth
        shmin_curve[depth_array < self.sf_depth_msl] = hydrostatic_pressure_curve[depth_array < self.sf_depth_msl] 

        return shmin_curve

    def add_scenario(self, scenario_name: str, **kwargs):

        if 'fluid_type' in kwargs and kwargs['fluid_type'] != self.fluid_type:
            # If fluid_type is provided, load the PVT data for the new fluid type
            pvt_data = load_pvt_data(self.pvt_path, kwargs['fluid_type'], load_fluid=True)
            kwargs['pvt_data'] = pvt_data
            kwargs['fluid_composition'] = pvt_data[kwargs['fluid_type']]['metadata']['composition']

            temperature_vector = pvt_data['temperature']
            pressure_vector = pvt_data['pressure']
            rho_matrix = pvt_data[kwargs['fluid_type']]['rho']

            kwargs['fluid_interpolator'] = RectBivariateSpline(pressure_vector, temperature_vector, rho_matrix)

        elif 'specific_gravity' in kwargs and kwargs['specific_gravity'] is not None:
            # If specific_gravity is provided, set fluid_type to None and fluid_interpolator to None
            kwargs['fluid_type'] = None
            kwargs['fluid_interpolator'] = None
            kwargs['pvt_data'] = None
            kwargs['fluid_composition'] = None
        # # Check that either fluid_type or specific_gravity is provided, but not both
        if 'fluid_type' not in kwargs and 'specific_gravity' not in kwargs:
            if self.fluid_type is None and self.specific_gravity is None:
                raise ValueError("Either fluid_type or specific_gravity should be provided, not both.")
            else:
                kwargs['fluid_type'] = self.fluid_type
                kwargs['specific_gravity'] = self.specific_gravity

        


        defaults = {
            'z_fluid_contact': self.z_fluid_contact,
            'p_fluid_contact': self.p_fluid_contact,
            'init_curves': self.init_curves,
            'brine_interpolator': self.brine_interpolator,
            'fluid_type': self.fluid_type,
            'fluid_composition': self.fluid_composition,
            'pvt_data': self.pvt_data,
            'fluid_interpolator': self.fluid_interpolator,
            'specific_gravity': self.specific_gravity,
            'z_resrv': self.z_resrv,
            'p_resrv': self.p_resrv,
            # Add other default values as needed
        }

        # Update defaults with user-provided values
        defaults.update(kwargs)



        # Store the scenario
        scenario = self.scenario_manager.create_scenario(name=scenario_name, **defaults)
        scenario.compute_pressure_profile()


    def manage_scenarios(self):
        # Check if the first scenario is 'None' and should be computed as hydrostatic
        # Skip the first entry in the input_scenarios dictionary

        if self.input_scenarios:
            scenarios_iter = iter(self.input_scenarios.items())
            next(scenarios_iter)  
            for scenario_name, scenario_pressure in scenarios_iter:
                print(f'{scenario_pressure=}')
                try:
                    p_delta = float(scenario_pressure)
                except:
                    p_delta = None  # Handle cases where pressure is not a float
                print(f'{scenario_name=} {scenario_pressure=}')
                self.add_scenario(scenario_name=scenario_name, p_delta=p_delta, from_resrvr=True,)



        
        elif self.z_fluid_contact is not None and self.default_hs_scenario:
            # Handle default hydrostatic scenario if no input scenarios are provided
            self.add_scenario(
                scenario_name='hydrostatic',
                from_resrvr=True,
            )
            

    
    def compute_barrier_leakage(self, well: Well, barrier_name: str) -> pd.DataFrame:
        """ 
        Compute leakage rate from the given barrier

        Args:
            well (Well): well information
            barrier_name (str): barrier to check the leakage rate
        Returns:
            pd.DataFrame: DataFrame containing leakage rates for different scenarios and permeabilities.

        """
        if well.inventory['barriers']:
            # for convenience
            barrier_perm = well.barrier_perm

            # Estimate CO2 leakage in [tons/day] after a trancient period
            sc_names = self.scenario_manager.scenarios.keys()

            # Initialize a DataFrame to store the leakage rates
            df = pd.DataFrame(columns=["p_brine_above_barrier", "p_fluid_below_barrier", "rho_brine_below_barrier", "rho_fluid_below_barrier"], index=sc_names)

            # Retrieve common init curves
            depth = self.init_curves['depth'].values
            temperature = self.init_curves['temperature'].values

            # barrier geometries
            barrier_props = well.compute_barrier_props(barrier_name)

            b_top   = barrier_props['top']
            b_bottom= barrier_props['bottom']
            
            # Retrieve temperature at the top and bottom of the barrier
            b_top_temp = np.interp(b_top, depth, temperature)
            b_bottom_temp = np.interp(b_bottom, depth, temperature)

            for sc_name in sc_names:
                # Retrieve and interpolate pressure and density values
                p_brine_ab, p_fluid_bb, rho_brine_ab, rho_fluid_bb = self._retrieve_and_interpolate_values(sc_name = sc_name, 
                                                                                                           top = b_top,
                                                                                                           bottom = b_bottom, 
                                                                                                           top_temperature = b_top_temp, 
                                                                                                           bottom_temperature = b_bottom_temp, 
                                                                                                           depth = depth)

                # Store retrieved values in the DataFrame
                df.loc[sc_name, "p_brine_above_barrier"] = p_brine_ab
                df.loc[sc_name, "p_fluid_below_barrier"] = p_fluid_bb
                df.loc[sc_name, "rho_brine_below_barrier"] = rho_brine_ab
                df.loc[sc_name, "rho_fluid_below_barrier"] = rho_fluid_bb

                # Check if the barrier has permeabilities
                try:
                    perms = barrier_perm['kv'].values()
                except Exception:
                    perms = barrier_perm['kv']

                # Compute leakage rates for different permeabilities and store in df
                for k in perms:
                    df[k] = np.nan

                    for idx, row in df.iterrows():
                        df.loc[idx, k] = leakage_proxy(rho_fluid_below_barrier = row['rho_fluid_below_barrier'],
                                                    rho_brine_below_barrier = row['rho_brine_below_barrier'],
                                                    p_fluid_below_barrier = row['p_fluid_below_barrier'],
                                                    p_brine_above_barrier = row['p_brine_above_barrier'],
                                                    permeability = k,
                                                    barrier_props = barrier_props)

            return df
        
        else:
            print(f'No barriers declared in well {well.header["well_name"]}')

    def _retrieve_and_interpolate_values(self, sc_name: str, top: float, bottom: float, top_temperature: float, bottom_temperature: float, depth: np.ndarray) -> tuple:
        """
        Retrieve and interpolate pressure and density values at the top and bottom of the barrier.

        Args:
            sc_name (str): Scenario name
            top (float): Top depth of the barrier
            bottom (float): Bottom depth of the barrier
            top_temperature (float): Temperature at the top of the barrier
            bottom_temperature (float): Temperature at the bottom of the barrier
            depth (np.ndarray): Depth array

        Returns:
            tuple: Interpolated pressure and density values (p_brine_above_barrier, p_fluid_below_barrier, rho_brine_above_barrier, rho_fluid_below_barrier)
        """
        fluid_pressure = self.scenario_manager.scenarios[sc_name].init_curves['fluid_pressure']
        hydrst_pressure = self.scenario_manager.scenarios[sc_name].init_curves['hydrostatic_pressure']
        fluid_interp = self.scenario_manager.scenarios[sc_name].fluid_interpolator
        brine_interp = self.scenario_manager.scenarios[sc_name].brine_interpolator

        p_fluid_below_barrier = np.interp(bottom, depth, fluid_pressure)
        p_brine_above_barrier = np.interp(top, depth, hydrst_pressure)

        rho_fluid_below_barrier = get_rho_from_pvt_data(pressure=p_fluid_below_barrier, temperature=bottom_temperature, rho_interpolator=fluid_interp)
        rho_brine_above_barrier = get_rho_from_pvt_data(pressure=p_brine_above_barrier, temperature=top_temperature, rho_interpolator=brine_interp)

        return p_brine_above_barrier, p_fluid_below_barrier, rho_brine_above_barrier, rho_fluid_below_barrier

    def scenarios_summary(self):
        return self.scenario_manager.get_scenarios_summary()




