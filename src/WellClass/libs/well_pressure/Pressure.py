from dataclasses import dataclass, field
from typing import Union, Dict
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d, RectBivariateSpline


from .PressureScenarioManager import PressureScenarioManager
from .helper_func import load_pvt_data, compute_hydrostatic_pressure

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

        brine_rho_matrix = self.pvt_data['brine']['rho']
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
            shmin_curve = np.copy(hydrostatic_pressure_curve) # Copy hydrostatic pressure as a starting point
            below_sf_mask = depth_array >= self.sf_depth_msl # Mask for depths below seafloor
            shmin_curve[below_sf_mask] += (depth_array[below_sf_mask] - self.sf_depth_msl) * SHMIN_FAC

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




        #  # Assign fluid_type from the class if it's not overridden in the method call

        # # Check that either fluid_type or specific_gravity is provided, but not both
        # if fluid_type is not None and specific_gravity is not None:
        #     raise ValueError("Either fluid_type or specific_gravity should be provided, not both.")



        # print(f'{scenario_name=}')
        # fluid_type = fluid_type or self.fluid_type
        # pvt_arrays = list(self.pvt_data.keys())



                
        

        # # Determine if specific gravity is provided for this scenario
        # specific_gravity_provided = specific_gravity is not None

        # # Ensure init_curves and z_fluid_contact are included in kwargs or add them if not present
        # kwargs.setdefault('init_curves', self.init_curves)
        # kwargs.setdefault('z_fluid_contact', self.z_fluid_contact)
        # kwargs.setdefault('fluid_type', fluid_type)
        # print(f'{specific_gravity=}')
        # print(f'{specific_gravity is not None=}')


        # # Handle the specific gravity and PVT data based on the scenario
        # if specific_gravity_provided:
        #     # Scenario added with specific gravity
        #     kwargs['specific_gravity'] = specific_gravity
        #     kwargs['pvt_data'] = None
        #     kwargs['fluid_interpolator'] = None
        # else:
        #     # Scenario added with fluid PVT data
        #     # Check if the fluid_type is different from the one already stored or if PVT data is not in kwargs
        #     if fluid_type != self.fluid_type or 'pvt_data' not in kwargs:
        #         # Load PVT data for the new fluid type
        #         pvt_data = load_pvt_data(self.pvt_path, fluid_type, load_fluid=True)
        #         kwargs['pvt_data'] = pvt_data
        #         temperature_vector = pvt_data['temperature']
        #         pressure_vector = pvt_data['pressure']
        #         rho_matrix = pvt_data[fluid_type]['rho']
        #         kwargs['fluid_interpolator'] = RectBivariateSpline(pressure_vector, temperature_vector, rho_matrix)
        #     else:
        #         # Use the existing PVT data and interpolators
        #         kwargs.setdefault('pvt_data', self.pvt_data)
        #         kwargs.setdefault('fluid_interpolator', self.fluid_interpolator)

        # # Always use the brine interpolator from the Pressure instance
        # kwargs.setdefault('brine_interpolator', self.brine_interpolator)

       
        # # Create a new PressureScenario instance with the given name and parameters
        # scenario = self.scenario_manager.create_scenario(name=scenario_name, **kwargs)
        # scenario.compute_pressure_profile()

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

    
