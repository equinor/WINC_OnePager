from collections import defaultdict

from dataclasses import dataclass
from pathlib import Path
import json

from typing import Union
import numpy as np
import pandas as pd
import math

from ..pvt.pvt import ( get_hydrostatic_P,
                       get_pvt)
from ..well_class.well_class import Well

from .barrier_pressure import (
    compute_barrier_leakage
)

from  .pressure_scenario import FluidP_scenario

from .co2_pressure import (
    _get_shmin,
    _integrate_pressure,
    _get_rho_in_pressure_column
)

from ..utils.compute_intersection import compute_intersection



from scipy.interpolate import RectBivariateSpline


@dataclass              # @dataclass(kw_only=True)
class Pressure:
    """ This is used to compute leakeage rate of legacy well

        Args:
            header (dict): well header
            reservoir_P (dict): reservoir pressure
            co2_datum (float): co2 datum depth
            pvt_path (str): directory where PVT files are located.            
    """
    header          : dict = None
    reservoir_P     : dict = None    
    co2_datum       : dict = None
    pvt_path        : Union[str, Path] = None
    barriers        : dict = None
    max_pressure_pos: Union[dict, list, float, int] = None
    pressure_scenarios : dict = None
    mixture_name   :  str = None
    mixture_composition   :  str = None
    default_hs_scenario : bool = True

    def __post_init__(self):
        self._get_mixture_info()
        self.init_pressure_curves()
        self.check_init_pressure()
        self._init_CO2_pressures_table()
        self.check_scenarios_from_input()


    # TODO(hzh): non-pure function!!!
    def _get_mixture_info(self):
        
        if isinstance(self.pvt_path, str):
            pvt_path = Path(self.pvt_path)
        else:
            pvt_path = self.pvt_path

        with open(pvt_path / "metadata.json", "r") as file:
            mixture_info = json.load(file)
        
        self.mixture_name = mixture_info['name']
        self.mixture_composition = mixture_info['composition']
        print(f'Computing pressures for {self.mixture_name} ({self.mixture_composition})')
        

    def init_pressure_curves(self):
        '''
        Initiates Pressure table for the wellbore
        By default creates a hydrostatic gradient curve and a Sh min curve
        '''

        #Calculate depth, temp(depth) hydrostatic_pressure(depth), H2ORHO(depth for hydrostatic pressure)
        init_curves = get_hydrostatic_P(self.header, pvt_path=self.pvt_path) 

        #Include Minimum horizontal stress
        init_curves = _get_shmin(self.header, init_curves)

        #Store tables in main pressure_CO2 table 
        self.init_curves = init_curves

    def check_init_pressure(self):
        '''
        Calculates hydrostatic pressure at reservoir_P['depth_msl'] and stores it in reservoir_P['hydrostatic_pressure']

        reservoir_pressure		
        depth_msl	RP1	RP2
        2238		        20

        then P_init['depth_msl'] is the reference depth 2238

        Also delta-pressures RP1 and RP2 are read and used later for leakage calculations. (self.reservoir_P[])

        Note: The interpretation of the numbers set for RP1 and RP2 is a bit unclear. 
              RP means reservoir pressure - but for the current implementation 
              it is interpreted as delta pressure: delta wrt hydrostatic pressure. 
              Initial implementation had the possibility to have absolute pressure AND a delta-pressure. But the delta-pressure had to be 
              given as a string "+ 20" or "- 15". Bit the + and - tended to create problems. 
              So alternatively one could distingwuish between RP and DP  if it is important to give in absolute pressure as well.
              But the most interesting input IS actually the change in pressure compared with hydrostatic pressure.

        '''
        # If no reservoir pressure is given, then we assume the reservoir is at the CO2 datum
        if self.reservoir_P is None:
            self.reservoir_P = {'depth_msl': self.co2_datum}

        #Get initial pressure table from reservoir_P dictionary
        P_init = self.reservoir_P

        # Get reference depth from reservoir_P, e.g. top reservoir
        ref_z        = P_init['depth_msl']                       

        
        #Hydrostatic pressure at reference depth (ref_z)
        ref_p        = np.interp(ref_z, self.init_curves['depth_msl'], self.init_curves[ 'hs_p'])
        print(f"Hydrostatic pressure at reference depth {ref_z:.0f} is {ref_p:.2f}")

        #Store value in reservoir_P table
        self.reservoir_P['hydrostatic_pressure'] = float(ref_p)



    def _init_pressure_scenarios_dict(self):
        '''
        Returns a dictionary with default values for a new pressure scenario
        '''
        # Define the default dictionary
        default_dict = {
            'name': None,
            'z_MSAD': None,
            'from_resrvr': None,
            'p_MSAD': None,
            'p_resrv': None,
            'z_resrv': None,
            'p_delta': None
        }
            
        self.pressure_scenarios = defaultdict(lambda: default_dict.copy())

    def _init_CO2_pressures_table(self):
         #Define and assign a multiindex that groups columns generic for all scenarios: depth_msl, temp, hs_p, RHOH2O and Shmin 
        index_init = pd.MultiIndex.from_product([['init'], self.init_curves.columns])
        self.pressure_CO2 = self.init_curves.copy()
        self.pressure_CO2.columns = index_init
        

    def check_scenarios_from_input(self):
        '''
        Reads input scenarios and creates a dictionary of scenarios to compute pressures for
        '''
        #Initialize dictionary of scenarios
        self._init_pressure_scenarios_dict()

        MAX_PRESSURE_NAME = 'max_p'

        #Get hydrostatic pressure at reservoir depth
        ref_p = self.reservoir_P['hydrostatic_pressure']

        #Compute number of scenarios
        n_scenarios = max(len(self.reservoir_P)-2, 0)

        #Initialize scenario counter
        scenario_counter = 0

        #If no scenarios are given, then we have at least one scenario: hydrostatic pressure
        if self.default_hs_scenario:
            self._add_hydrostatic_scenario(scenario_counter, ref_p)

            #Increment scenario counter
            scenario_counter += 1
        
        #Iterate over input scenarios
        if n_scenarios > 0:
            scenario_counter = self._add_input_scenarios(scenario_counter, ref_p)

        #Iterate over maximum pressure scenarios
        if not self.max_pressure_pos is None:
            self._add_max_pressure_scenarios(scenario_counter, MAX_PRESSURE_NAME)



    def _add_hydrostatic_scenario(self, scenario_counter, ref_p):
        '''
        Include hydrostatic pressure scenario to pressure_scenarios dictionary
        '''
        self.pressure_scenarios[scenario_counter]['name'] = 'hydrostatic'
        self.pressure_scenarios[scenario_counter]['p_resrv'] = ref_p
        self.pressure_scenarios[scenario_counter]['from_resrvr'] = True
        self.pressure_scenarios[scenario_counter]['z_resrv'] = self.reservoir_P['depth_msl']
        self.pressure_scenarios[scenario_counter]['p_delta'] = 0
        self.pressure_scenarios[scenario_counter]['z_co2_datum'] = self.co2_datum

        sc_pressure = self._compute_scenario_profiles(self.pressure_scenarios[scenario_counter])

        self.pressure_scenarios[scenario_counter]['p_MSAD'] = sc_pressure.p_MSAD
        self.pressure_scenarios[scenario_counter]['z_MSAD'] = sc_pressure.z_MSAD

    def _add_input_scenarios(self, scenario_counter, ref_p):
        '''
        If more scenarios are included in input data, then they are added to the pressure_scenarios dictionary
        '''
        keys_to_skip = {'depth_msl', 'hydrostatic_pressure'}
        first_as_hydrostatic = False

        for sc_name, sc_pressure in self.reservoir_P.items():
            if sc_name in keys_to_skip:
                continue

            if sc_pressure is None or (isinstance(sc_pressure, float) and math.isnan(sc_pressure)):
                if not first_as_hydrostatic:

                    self.pressure_scenarios[0]['name'] = sc_name


                    self.pressure_CO2.columns = self.pressure_CO2.columns.set_levels(
                        self.pressure_CO2.columns.levels[0].str.replace('hydrostatic', sc_name), level=0)
                    
                    first_as_hydrostatic = True
                continue

            magnitude = self._parse_pressure_magnitude(sc_pressure)
            p_resrv = ref_p + magnitude
            self.pressure_scenarios[scenario_counter]['name'] = sc_name
            self.pressure_scenarios[scenario_counter]['p_resrv'] = p_resrv
            self.pressure_scenarios[scenario_counter]['z_resrv'] = self.reservoir_P['depth_msl']
            self.pressure_scenarios[scenario_counter]['from_resrvr'] = True
            self.pressure_scenarios[scenario_counter]['p_delta'] = magnitude

            sc_pressure = self._compute_scenario_profiles(self.pressure_scenarios[scenario_counter])

            self.pressure_scenarios[scenario_counter]['p_MSAD'] = sc_pressure.p_MSAD
            self.pressure_scenarios[scenario_counter]['z_MSAD'] = sc_pressure.z_MSAD



            scenario_counter += 1

        return scenario_counter

    def _parse_pressure_magnitude(self, sc_pressure):
        '''
        Parses pressure magnitude from string to float
        '''
        try:
            if isinstance(sc_pressure, str):
                sc_pressure = sc_pressure.replace(" ", "")
            return float(sc_pressure)
        except Exception:
            return sc_pressure
    
    def _add_max_pressure_scenarios(self, scenario_counter, MAX_PRESSURE_NAME):
        '''
        Includes maximum pressure scenarios in the pressure_scenarios dictionary
        Maximum pressure scenarios are calculated from the Shmin values at specific depths above the reservoir
        '''
        # Determine the iterable based on the type of self.max_pressure_pos
        # Check if max_pressure_pos is a dictionary of barriers
        if isinstance(self.max_pressure_pos, dict):  
            print(f'max_pressure_pos is a dictionary of barriers')
            iterable = [(self.max_pressure_pos['bottom_msl'][idx], key) for idx, key in self.max_pressure_pos['barrier_name'].items()]

        # Check if max_pressure_pos is a depth value or a list of depths
        elif isinstance(self.max_pressure_pos, (list, float, int)):
            print(f'max_pressure_pos is a value')
            if isinstance(self.max_pressure_pos, (float, int)):  # Make it a list with one element
                self.max_pressure_pos = [self.max_pressure_pos]
            iterable = [(depth, f"at_{int(depth)}") for depth in self.max_pressure_pos]
        else:
            raise ValueError("Invalid type for max_pressure_pos")

        # Process the iterable
        for barr_depth, key in iterable:
            sc_name = f"{MAX_PRESSURE_NAME}_{key}"
            self.pressure_scenarios[scenario_counter]['name'] = sc_name
            self.pressure_scenarios[scenario_counter]['z_MSAD'] = barr_depth
            self.pressure_scenarios[scenario_counter]['from_resrvr'] = False

            sc_pressure = self._compute_scenario_profiles(self.pressure_scenarios[scenario_counter])

            self.pressure_scenarios[scenario_counter]['p_MSAD'] = sc_pressure.p_MSAD
            self.pressure_scenarios[scenario_counter]['p_resrv'] = sc_pressure.p_resrv
            self.pressure_scenarios[scenario_counter]['z_resrv'] = sc_pressure.z_resrv 
            self.pressure_scenarios[scenario_counter]['p_delta'] = sc_pressure.p_delta

            scenario_counter += 1


    def _compute_scenario_profiles(self, pressure_scenario: dict):
        #Retrieve pressure, temperature and density fields for CO2 and H2O
        pvt_T, pvt_P, pvt_RHO_CO2, pvt_RHO_H2O =  get_pvt(self.pvt_path)        

        #An interpolator. Used later to look-up densities given pressure p and temperature t
        get_rho_h2o = RectBivariateSpline(pvt_P, pvt_T, pvt_RHO_H2O)
        get_rho_co2 = RectBivariateSpline(pvt_P, pvt_T, pvt_RHO_CO2)
        
        sc_pressure = FluidP_scenario(**pressure_scenario)
        sc_pressure.compute_pressure_profile(init_table = self.init_curves,
                                             well_header = self.header,
                                             rho_co2_getter = get_rho_co2,
                                             rho_h2o_getter = get_rho_h2o)
        
        

        self._merge_scenario_profiles(sc_pressure)

        return sc_pressure
        



    def _merge_scenario_profiles(self, sc_pressure: FluidP_scenario):

        #Drop repeated columns
        sc_pressure.P_table = sc_pressure.P_table.drop(columns = self.pressure_CO2['init'].columns)

        #Setup multiindex
        index_table = pd.MultiIndex.from_product([[sc_pressure.name], sc_pressure.P_table.columns])
        sc_pressure.P_table.columns = index_table
            
        #Concatenate init table with scenario table
        self.pressure_CO2 = pd.concat([self.pressure_CO2, sc_pressure.P_table], axis= 1)




    def compute_barrier_leakage(self, well: Well, barrier_name: str) -> pd.DataFrame:
        """ Compute leakage rate from the given barrier

            Args:
                well (Well): well information
                barrier_name (str): barrier to check the leakage rate
        """
        if well.inventory['barriers']:
            # for convenience
            barrier_perm = well.barrier_perm

            # barrier geometries
            barrier_props = well.compute_barrier_props(barrier_name)

            # Estimate CO2 leakage in [tons/day] after a trancient period
            barrier_leakage = compute_barrier_leakage(barrier_perm, self.pressure_scenarios, self.pressure_CO2, barrier_props)

            return barrier_leakage
        
        else:
            print(f'No barriers declared in well {well.header["well_name"]}')

    def create_pressure_scenario(self, name: str = None, z_MSAD: float = None, from_resrvr: bool = None, p_MSAD: float = None, p_resrv: float = None, z_resrv: float = None, p_delta: float = None):
        scenario_counter = len(self.pressure_scenarios)
        self.pressure_scenarios[scenario_counter] = {
            'name': name,
            'z_MSAD': z_MSAD,
            'from_resrvr': from_resrvr,
            'p_MSAD': p_MSAD,
            'p_resrv': p_resrv,
            'z_resrv': z_resrv,
            'p_delta': p_delta,
            'z_co2_datum': self.co2_datum
        }

        if name is None:
            # Return error ir scenario name is not provided
            raise ValueError("The 'name' parameter is required.")
        
        if from_resrvr is None:
            # Return error if 'from_resrvr' is not provided
            raise ValueError("The 'from_resrvr' parameter is required.")
        
        if from_resrvr:
            # Compute pressure depending on the input parameters provided
            if p_delta is None and (p_resrv is None or z_resrv is None):
                # Return error if 'p_delta' or both 'p_resrv' and 'z_resrv' are not provided
                raise ValueError("If 'from_resrvr' is True, you must provide either 'p_delta' or both 'p_resrv' and 'z_resrv'.")
            
            if p_resrv is not None and z_resrv is not None:
                # Compute pressure if both 'p_resrv' and 'z_resrv' are provided
                # Interpolate hydrostatic pressure at z_resrv
                hydrostatic_pressure = np.interp(z_resrv, self.init_curves['depth_msl'], self.init_curves['hs_p'])

                # Update z_co2_datum with provided reservoir depth
                self.pressure_scenarios[scenario_counter]['z_co2_datum'] = z_resrv

                
                # Compute reservoir pressure
                self.pressure_scenarios[scenario_counter]['p_delta'] = p_resrv - hydrostatic_pressure

            elif p_delta is not None and z_resrv is not None:
                # Compute pressure if 'p_delta' and 'z_resrv' are provided
                                
                # Interpolate hydrostatic pressure at z_resrv
                hydrostatic_pressure = np.interp(z_resrv, self.init_curves['depth_msl'], self.init_curves['hs_p'])

                # Compute reservoir pressure
                self.pressure_scenarios[scenario_counter]['p_resrv'] = hydrostatic_pressure + p_delta

            elif p_delta is not None:
                # Compute pressure if only 'p_delta' is provided.
                # Use self.co2_datum as z_resrv
                
                z_resrv = self.co2_datum
                hydrostatic_pressure = np.interp(z_resrv, self.init_curves['depth_msl'], self.init_curves['hs_p'])
                self.pressure_scenarios[scenario_counter]['p_resrv'] = hydrostatic_pressure + p_delta
                self.pressure_scenarios[scenario_counter]['z_resrv'] = z_resrv
            else:
                raise ValueError("Invalid combination of parameters for 'from_resrvr' = True.")
        else:
            if z_MSAD is None:
                raise ValueError("If 'from_resrvr' is False, you must provide 'z_MSAD'.")


        sc_pressure = self._compute_scenario_profiles(self.pressure_scenarios[scenario_counter])

        self.pressure_scenarios[scenario_counter]['p_MSAD'] = sc_pressure.p_MSAD
        self.pressure_scenarios[scenario_counter]['z_MSAD'] = sc_pressure.z_MSAD
        self.pressure_scenarios[scenario_counter]['p_delta'] = sc_pressure.p_delta
        self.pressure_scenarios[scenario_counter]['p_resrv'] = sc_pressure.p_resrv
        self.pressure_scenarios[scenario_counter]['z_resrv'] = sc_pressure.z_resrv
        self.pressure_scenarios[scenario_counter]['z_co2_datum'] = sc_pressure.z_co2_datum


    @property
    def to_json(self):
        return json.dumps(self.__dict__, indent=4)