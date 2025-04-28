from dataclasses import dataclass, field
from typing import Union, Dict
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
        collated_profiles (pd.DataFrame, optional): Collated profiles for all scenarios.

    """


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
    fluid_composition: str = None
    specific_gravity: float = None
    pvt_data: Dict[str, Dict[str, np.ndarray]] = field(init=False)
    brine_interpolator: RectBivariateSpline = field(init=False)
    fluid_interpolator: RectBivariateSpline = field(init=False)
    ip_shmin_data: np.ndarray = field(default=None) # User-provided Shmin data
    init_curves: pd.DataFrame = field(init=False)
    scenario_manager: PressureScenarioManager = field(init=False)
    input_scenarios: Dict[str, Union[float, str, None]] = field(default_factory=dict)
    default_hs_scenario: bool = True
    collated_profiles: pd.DataFrame = field(init=None)
    salinity: float = 3.5  # Salinity in percentage (default is 3.5% for seawater)
    shmin_gradient: float = SHMIN_FAC  # Gradient for Shmin calculation


    def __post_init__(self):
        """
        Initializes PVT data, interpolators, initial curves, and manages scenarios 
        after class instantiation.
        """

        if self.specific_gravity is not None:
            # Load brine data only
            self.pvt_data = load_pvt_data(self.pvt_path, load_fluid=False)
            # Initialize brine interpolator only
            self._initialize_interpolators(load_fluid=False)
            self.fluid_composition = None
        else:
            # Load both brine and fluid data
            self.pvt_data = load_pvt_data(self.pvt_path, self.fluid_type, load_fluid=True)
            self.fluid_composition = self.pvt_data[self.fluid_type]['metadata']['composition']
            # Initialize both brine and fluid interpolators
            self._initialize_interpolators(load_fluid=True)


        self.init_curves = self._compute_init_curves()
        self.scenario_manager = PressureScenarioManager()
        self.manage_scenarios()

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
        min_horizontal_stress[min_horizontal_stress<hydrostatic_pressure_curve] = hydrostatic_pressure_curve[min_horizontal_stress<hydrostatic_pressure_curve]  
        # min_horizontal_stress[hydrostatic_pressure_curve>min_horizontal_stress] = hydrostatic_pressure_curve[hydrostatic_pressure_curve>min_horizontal_stress]
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
        if load_fluid:
            fluid_rho_matrix = self.pvt_data[self.fluid_type]['rho']
            self.fluid_interpolator = RectBivariateSpline(pressure_vector, temperature_vector, fluid_rho_matrix)


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
            depth_values, shmin_values = self.ip_shmin_data.T
            shmin_interpolator = interp1d(depth_values, shmin_values, bounds_error=False, fill_value="extrapolate")

            shmin_curve = shmin_interpolator(depth_array)

        else:
            # No user data provided, use the empirical formula
            
            # interpolate hydrostatic pressure at mudline depth (seafloor depth)
            pressure_ml = np.interp(self.sf_depth_msl, depth_array,  hydrostatic_pressure_curve)
            
            depth_ml = depth_array - self.sf_depth_msl  # depth below mean sea level

            shmin_curve = pressure_ml + depth_ml*self.shmin_gradient
            shmin_curve[shmin_curve<0] = hydrostatic_pressure_curve[shmin_curve<0]


        return shmin_curve

    def add_scenario(self, scenario_name: str, fluid_type: str = None, specific_gravity: float = None, **kwargs):
        # Check that either fluid_type or specific_gravity is provided, but not both
        if fluid_type is None and specific_gravity is None:
            fluid_type = self.fluid_type
            specific_gravity = self.specific_gravity
            if fluid_type is not None and specific_gravity is not None:
                raise ValueError("Either fluid_type or specific_gravity should be provided, not both.")        


        # If fluid_type is provided, use it and ensure specific_gravity is not used
        if fluid_type is not None:
            # Ensure specific_gravity is not set in kwargs
            kwargs['specific_gravity'] = None

            if fluid_type != self.fluid_type or 'pvt_data' not in kwargs:
                # Load PVT data for the new fluid type
                pvt_data = load_pvt_data(self.pvt_path, fluid_type, load_fluid=True)
                kwargs['pvt_data'] = pvt_data
                kwargs['fluid_type'] = fluid_type
                kwargs['fluid_composition'] = pvt_data[fluid_type]['metadata']['composition']


                temperature_vector = pvt_data['temperature']
                pressure_vector = pvt_data['pressure']
                rho_matrix = pvt_data[fluid_type]['rho']

                kwargs['fluid_interpolator'] = RectBivariateSpline(pressure_vector, temperature_vector, rho_matrix)
                

            else:
                # Use the existing PVT data and interpolators
                kwargs.setdefault('pvt_data', self.pvt_data)
                kwargs.setdefault('fluid_interpolator', self.fluid_interpolator)
                kwargs.setdefault('fluid_composition', self.fluid_composition)
                kwargs.setdefault('fluid_type', self.fluid_type)
            
        elif specific_gravity is not None:
            kwargs['specific_gravity'] = specific_gravity
            kwargs['pvt_data'] = None
            kwargs['fluid_interpolator'] = None
            # Ensure fluid_type is not set in kwargs
            kwargs['fluid_type'] = None




        # Ensure init_curves and z_fluid_contact are included in kwargs or add them if not present
        kwargs.setdefault('init_curves', self.init_curves)
        # kwargs.setdefault('z_fluid_contact', self.z_fluid_contact)

        # Always use the brine interpolator from the Pressure instance
        kwargs.setdefault('brine_interpolator', self.brine_interpolator)

        # Create the PressureScenario with the parameters and store it
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
                    fluid_type=self.fluid_type,
                    z_fluid_contact=self.z_fluid_contact,
                    p_delta=0  # Assuming hydrostatic pressure
                )

            
            # Parse the remaining scenarios and create PressureScenario instances
            for sc_name, sc_pressure in scenarios_iter:
                if sc_pressure is not None and isinstance(sc_pressure, str):
                    # Convert string to float, handling any special string cases here

                    p_delta = float(sc_pressure)
                    print(f'{p_delta=}')
                    self.add_scenario(
                        scenario_name=sc_name,
                        from_resrvr=True,
                        fluid_type=self.fluid_type,
                        z_fluid_contact=self.z_fluid_contact,
                        p_delta=p_delta
                    )

        
        elif self.z_fluid_contact is not None and self.default_hs_scenario:
            # Handle default hydrostatic scenario if no input scenarios are provided
            self.add_scenario(
                scenario_name='hydrostatic',
                from_resrvr=True,
                z_fluid_contact=self.z_fluid_contact,
                fluid_type=self.fluid_type,
            )
            
    def collate_all_profiles(self):
        # Call the collate_scenario_profiles method from the scenario manager
        common_data = self.init_curves[['depth', 'temperature', 'hydrostatic_pressure', 'min_horizontal_stress']]
        self.collated_profiles = self.scenario_manager.collate_scenario_profiles(common_data)

    
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


