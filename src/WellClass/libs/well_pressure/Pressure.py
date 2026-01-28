import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, Union

import numpy as np
import pandas as pd
import scipy.constants as const
from scipy.interpolate import RectBivariateSpline, interp1d

from ..pvt.pvt import get_rho_from_pvt_data
from ..well_class.well_class import Well
from .barrier_pressure import leakage_proxy
from .PressureScenarioManager import PressureScenarioManager
from .PVTDataManager import PVTDataManager
from .WellConditions import WellConditions

# Constants
SHMIN_NAME = "Shmin"
SHMIN_FAC = 0.1695  # empirical factor for Shmin calculation


@dataclass
class Pressure:
    """
    Manages pressure scenarios for a legacy well based on geothermal gradient and fluid properties.

    This class orchestrates the computation of pressure profiles for various scenarios,
    using WellConditions for initial curves and PVTDataManager for fluid properties.

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
        specific_gravity (float, optional): Specific gravity if using constant density fluid.
        rho_brine (float, optional): Density of brine in kg/m^3, if known.
        ip_shmin_data (np.ndarray, optional): User-provided minimum horizontal stress data.
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
    rho_brine: Optional[float] = None  # Density of brine in kg/m^3, if known
    ip_shmin_data: np.ndarray = field(default=None)  # User-provided Shmin data
    input_scenarios: Dict[str, Optional[Union[float, str]]] = field(default_factory=dict)
    default_hs_scenario: bool = True
    salinity: float = 3.5  # Salinity in percentage (default is 3.5% for seawater)
    shmin_gradient: float = SHMIN_FAC  # Gradient for Shmin calculation

    # Computed attributes
    pvt_manager: PVTDataManager = field(init=False)
    well_conditions: WellConditions = field(init=False)
    pvt_data: Dict[str, Dict[str, np.ndarray]] = field(init=False)
    brine_interpolator: RectBivariateSpline = field(init=False)
    fluid_interpolator: RectBivariateSpline = field(init=False)
    init_curves: pd.DataFrame = field(init=False)
    scenario_manager: PressureScenarioManager = field(init=False)

    def __post_init__(self):
        """
        Initializes PVT data, interpolators, initial curves, and manages scenarios
        after class instantiation.
        """
        # Initialize PVT data manager
        self.pvt_manager = PVTDataManager(
            pvt_path=self.pvt_path,
            fluid_type=self.fluid_type,
            specific_gravity=self.specific_gravity,
            salinity=self.salinity,
        )

        # Set fluid composition from PVT manager
        self.fluid_composition = self.pvt_manager.fluid_composition

        # Get PVT data and interpolators
        self.pvt_data = self.pvt_manager.get_pvt_data()
        self.brine_interpolator = self.pvt_manager.get_brine_interpolator()
        self.fluid_interpolator = self.pvt_manager.get_fluid_interpolator()

        # Initialize well conditions
        self.well_conditions = WellConditions(
            sf_depth_msl=self.sf_depth_msl,
            well_td_rkb=self.well_td_rkb,
            well_rkb=self.well_rkb,
            sf_temp=self.sf_temp,
            geo_tgrad=self.geo_tgrad,
            pvt_data=self.pvt_data,
            brine_interpolator=self.brine_interpolator,
            rho_brine=self.rho_brine,
            salinity=self.salinity,
            shmin_gradient=self.shmin_gradient,
            ip_shmin_data=self.ip_shmin_data,
        )

        # Get initial curves from well conditions
        self.init_curves = self.well_conditions.get_curves()

        # Initialize scenario manager and scenarios
        self.scenario_manager = PressureScenarioManager()
        self.manage_scenarios()

    def add_scenario(self, scenario_name: str, **kwargs):
        """
        Add a new pressure scenario to the scenario manager.

        Args:
            scenario_name: Name of the scenario.
            **kwargs: Scenario-specific parameters.
        """
        # Handle different fluid type if provided
        if "fluid_type" in kwargs and kwargs["fluid_type"] != self.fluid_type:
            fluid_interpolator, fluid_composition, pvt_data = self.pvt_manager.load_additional_fluid(kwargs["fluid_type"])
            kwargs["pvt_data"] = pvt_data
            kwargs["fluid_composition"] = fluid_composition
            kwargs["fluid_interpolator"] = fluid_interpolator

        elif "specific_gravity" in kwargs and kwargs["specific_gravity"] is not None:
            # If specific_gravity is provided, set fluid_type to None and fluid_interpolator to None
            kwargs["fluid_type"] = None
            kwargs["fluid_interpolator"] = None
            kwargs["fluid_composition"] = None

        # Check that either fluid_type or specific_gravity is provided
        if "fluid_type" not in kwargs and "specific_gravity" not in kwargs:
            if self.fluid_type is None and self.specific_gravity is None:
                raise ValueError("Either fluid_type or specific_gravity should be provided, not both.")
            else:
                kwargs["fluid_type"] = self.fluid_type
                kwargs["specific_gravity"] = self.specific_gravity

        # Set defaults from class attributes
        defaults = {
            "z_fluid_contact": self.z_fluid_contact,
            "p_fluid_contact": self.p_fluid_contact,
            "init_curves": self.init_curves,
            "brine_interpolator": self.brine_interpolator,
            "fluid_type": self.fluid_type,
            "fluid_composition": self.fluid_composition,
            "pvt_data": self.pvt_data,
            "fluid_interpolator": self.fluid_interpolator,
            "specific_gravity": self.specific_gravity,
            "z_resrv": self.z_resrv,
            "p_resrv": self.p_resrv,
            "rho_brine": self.rho_brine,
        }

        # Update defaults with user-provided values
        defaults.update(kwargs)

        # Create and compute the scenario
        scenario = self.scenario_manager.create_scenario(name=scenario_name, **defaults)
        scenario.compute_pressure_profile()

    def manage_scenarios(self):
        # Check if the first scenario is 'None' and should be computed as hydrostatic
        # Skip the first entry in the input_scenarios dictionary

        if self.input_scenarios:
            scenarios_iter = iter(self.input_scenarios.items())
            next(scenarios_iter)
            for scenario_name, scenario_pressure in scenarios_iter:
                print(f"{scenario_pressure=}")
                try:
                    p_delta = float(scenario_pressure)
                except (ValueError, TypeError):
                    p_delta = None  # Handle cases where pressure is not a float
                print(f"{scenario_name=} {scenario_pressure=}")
                self.add_scenario(
                    scenario_name=scenario_name,
                    p_delta=p_delta,
                    from_resrvr=True,
                )

        elif self.z_fluid_contact is not None and self.default_hs_scenario:
            # Handle default hydrostatic scenario if no input scenarios are provided
            self.add_scenario(
                scenario_name="hydrostatic",
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
        if well.inventory["barriers"]:
            # for convenience
            barrier_perm = well.barrier_perm

            # Estimate CO2 leakage in [tons/day] after a trancient period
            sc_names = self.scenario_manager.scenarios.keys()

            # Initialize a DataFrame to store the leakage rates
            df = pd.DataFrame(
                columns=["p_brine_above_barrier", "p_fluid_below_barrier", "rho_brine_below_barrier", "rho_fluid_below_barrier"], index=sc_names
            )

            # Retrieve common init curves
            depth = self.init_curves["depth"].values
            temperature = self.init_curves["temperature"].values

            # barrier geometries
            barrier_props = well.compute_barrier_props(barrier_name)

            b_top = barrier_props["top"]
            b_bottom = barrier_props["bottom"]

            # Retrieve temperature at the top and bottom of the barrier
            b_top_temp = np.interp(b_top, depth, temperature)
            b_bottom_temp = np.interp(b_bottom, depth, temperature)

            for sc_name in sc_names:
                # Retrieve and interpolate pressure and density values
                p_brine_ab, p_fluid_bb, rho_brine_ab, rho_fluid_bb = self._retrieve_and_interpolate_values(
                    sc_name=sc_name, top=b_top, bottom=b_bottom, top_temperature=b_top_temp, bottom_temperature=b_bottom_temp
                )

                # Store retrieved values in the DataFrame
                df.loc[sc_name, "p_brine_above_barrier"] = p_brine_ab
                df.loc[sc_name, "p_fluid_below_barrier"] = p_fluid_bb
                df.loc[sc_name, "rho_brine_below_barrier"] = rho_brine_ab
                df.loc[sc_name, "rho_fluid_below_barrier"] = rho_fluid_bb

                # Check if the barrier has permeabilities
                try:
                    perms = barrier_perm["kv"].values()
                except Exception:
                    perms = barrier_perm["kv"]

                # Compute leakage rates for different permeabilities and store in df
                for k in perms:
                    df[k] = np.nan

                    for idx, row in df.iterrows():
                        df.loc[idx, k] = leakage_proxy(
                            rho_fluid_below_barrier=row["rho_fluid_below_barrier"],
                            rho_brine_below_barrier=row["rho_brine_below_barrier"],
                            p_fluid_below_barrier=row["p_fluid_below_barrier"],
                            p_brine_above_barrier=row["p_brine_above_barrier"],
                            permeability=k,
                            barrier_props=barrier_props,
                        )

            return df

        else:
            print(f"No barriers declared in well {well.header['well_name']}")

    def _retrieve_and_interpolate_values(
        self,
        sc_name: str,
        top: float,
        bottom: float,
        top_temperature: float,
        bottom_temperature: float,
    ) -> tuple:
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
        fluid_pressure = self.scenario_manager.scenarios[sc_name].init_curves["fluid_pressure"].to_numpy()
        hydrst_pressure = self.scenario_manager.scenarios[sc_name].init_curves["hydrostatic_pressure"].to_numpy()
        depth = self.scenario_manager.scenarios[sc_name].init_curves["depth"].to_numpy()
        fluid_interp = self.scenario_manager.scenarios[sc_name].fluid_interpolator
        brine_interp = self.scenario_manager.scenarios[sc_name].brine_interpolator

        p_fluid_below_barrier = np.interp(bottom, depth, fluid_pressure)
        p_brine_above_barrier = np.interp(top, depth, hydrst_pressure)

        rho_fluid_below_barrier = get_rho_from_pvt_data(pressure=p_fluid_below_barrier, temperature=bottom_temperature, rho_interpolator=fluid_interp)
        rho_brine_above_barrier = get_rho_from_pvt_data(pressure=p_brine_above_barrier, temperature=top_temperature, rho_interpolator=brine_interp)

        return p_brine_above_barrier, p_fluid_below_barrier, rho_brine_above_barrier, rho_fluid_below_barrier

    def scenarios_summary(self) -> pd.DataFrame:
        """
        Get a summary of all scenarios.

        Returns:
            pd.DataFrame: Summary table of all scenarios.
        """
        return self.scenario_manager.get_scenarios_summary()
