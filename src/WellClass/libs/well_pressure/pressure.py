
from dataclasses import dataclass
import json

from typing import Union
import numpy as np
import pandas as pd

from ..pvt.pvt import get_hydrostatic_P
from ..well_class.well_class import Well

from .barrier_pressure import (
    compute_barrier_leakage
)
from .co2_pressure import (
    compute_CO2_pressures,
    compute_MSAD
)

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
    pvt_path        : str  = None
    barriers        : dict = None
    max_pressure_pos: Union[dict, list, float, int] = None

    def __post_init__(self):
        self._check_init_pressure()
        self.pressure_CO2 = compute_CO2_pressures(self.header, self.reservoir_P, self.co2_datum, pvt_path=self.pvt_path, max_pressure_pos = self.max_pressure_pos)
        self.MSAD = compute_MSAD(p_init = self.reservoir_P, pt_df = self.pressure_CO2)

    # TODO(hzh): non-pure function!!!
    def _check_init_pressure(self):
        '''
        Calculates hydrostatic pressure at 
        reservoir_P is the entry looking something like this:

        reservoir_pressure		
        depth_msl	RP1	RP2
        2238		        20

        then P_init['depth_msl'] is the reference depth 2238

        Also delta-pressures RP1 and RP2 are read and used later for leakage calculations. (self.reservoir_P[])

        Note: The interpretation of the numbers set for RP1 and RP2 is a bit unclear. RP means reservoir pressure - but for the current
              implementation it is interpreted as delta pressure: delta wrt hydrostatic pressure. 
              Initial implementation had the possibility to have absolute pressure AND a delta-pressure. But the delta-pressure had to be 
              given as a string "+ 20" or "- 15". Bit the + and - tended to create problems. 
              So alternatively one could distingwuish between RP and DP  if it is important to give in absolute pressure as well.
              But the most interesting input IS actually the change in pressure compared with hydrostatic pressure.

        '''
        P_init = self.reservoir_P

        #Get hydrostatic pressure
        ref_z        = P_init['depth_msl']                                                  #A reference depth given in the input, e.g. top reservoir
        hydrostaticP = get_hydrostatic_P(self.header, pvt_path=self.pvt_path)               #Calculate depth, temp(depth) hydrostatic_pressure(depth), H2ORHO(depth for hydrostatic pressure)
        ref_p        = np.interp(ref_z, hydrostaticP['depth_msl'], hydrostaticP['hs_p'])    #Hydrostatic pressure at that depth
        print(f"Hydrostatic pressure at reference depth {ref_z:.0f} is {ref_p:.2f}")

        self.reservoir_P['hydrostatic_pressure'] = ref_p

        for key in P_init.copy():   #Should just loop the RPs and check for isnan for each instead of special treatment of RP1. Call it DP1 DP2 instead. Delta Pressure
            if key == 'RP1':
                RP1 = P_init[key]
                if RP1 is None or np.isnan(RP1):
                    print(f'RP1 set as hydrostatic P = {ref_p:.2f} bar')
                    self.reservoir_P['RP1'] = ref_p
                else:
                    print(f'RP1 is set as delta pressure, which yields P = {ref_p:.2f} {RP1:+.2f} = {ref_p + RP1:.2f} bar')
                    self.reservoir_P['RP1'] = ref_p + RP1
            elif key.startswith('RP'):
                RP = P_init[key]
                try:
                    if isinstance(RP, str):
                        RP = RP.replace(" ", "")
                    RP = float(RP)
                except Exception:
                    pass
                if isinstance(RP, float) or isinstance(RP, int):
                    # p = ref_p + RP                          #self.reservoir_P['RP1'] + RP
                    print(f'{key} is set as delta pressure, which yields P = {ref_p:.2f} {RP:+.2f} = {ref_p + RP:.2f} bar')
                    self.reservoir_P[key] = ref_p + RP
                else:
                    P_init.pop(key)
                    print(RP, 'ignored')
                    continue 



    def compute_barrier_leakage(self, well: Well, barrier_name: str) -> pd.DataFrame:
        """ Compute leakage rate from the given barrier

            Args:
                well (Well): well information
                barrier_name (str): barrier to check the leakage rate
        """

        # for convenience
        barrier_perm = well.barrier_perm

        # barrier geometries
        barrier_props = well.compute_barrier_props(barrier_name)

        # Estimate CO2 leakage in [tons/day] after a trancient period
        barrier_leakage = compute_barrier_leakage(barrier_perm, self.reservoir_P, self.pressure_CO2, barrier_props)

        return barrier_leakage

    @property
    def to_json(self):
        return json.dumps(self.__dict__, indent=4)

