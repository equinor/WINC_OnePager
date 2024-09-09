from dataclasses import dataclass
import json
from typing import Union
import numpy as np
import pandas as pd
from scipy.interpolate import RectBivariateSpline

from ..pvt.pvt import ( get_hydrostatic_P,
                       get_pvt)
from ..well_class.well_class import Well
from .barrier_pressure import (
    compute_barrier_leakage
)
from .co2_pressure import (
    _get_shmin,
    _integrate_pressure,
    _get_rho_in_pressure_column
)
from ..utils.compute_intersection import compute_intersection





@dataclass
class FluidP_scenario:
    name : str = None
    z_MSAD : float = np.nan
    from_resrvr : bool = None
    p_MSAD : float = np.nan
    p_resrv: float = np.nan 
    z_resrv: float = np.nan
    p_delta: float = np.nan


    def compute_pressure_profile(self,
                                 init_table: pd.DataFrame,
                                 well_header: dict,
                                 rho_co2_getter: RectBivariateSpline,
                                 rho_h2o_getter: RectBivariateSpline,
                                 z_co2_datum: float):
        """
        Depending on input defines which method to run.
        """
        if not isinstance(self.from_resrvr, bool):
            # class will compute max pressure from reservoir upwards
            raise ValueError("Invalid input parameters for computing pressure profile")
        elif self.from_resrvr:
            self.integrate_upwards(init_table, well_header,rho_co2_getter, rho_h2o_getter, z_co2_datum)
        else:
            self.integrate_downwards(init_table, well_header,rho_co2_getter, rho_h2o_getter, z_co2_datum)
            

    def integrate_upwards(self,
                          init_table: pd.DataFrame,
                          well_header: dict,
                          rho_co2_getter: RectBivariateSpline,
                          rho_h2o_getter: RectBivariateSpline,
                          z_co2_datum: float):
        """
        Method to compute pressure tables starting from reservoir pressure.
        CO2 pressure is integrated upwards from CO2 datum and input reservoir pressure (p_resrv).
        H2O pressure is integrated upwards and downwards from CO2 datum and input reservoir pressure (p_resrv).
        """
        print(f'Pressure scenario {self.name}: Compute pressures assuming reservoir pressure is {self.p_resrv:.2f} bar at {self.z_resrv:.2f} mTVDMSL')

        #Water
        p0 = self.p_resrv
        p_delta = self.p_delta


        water_p_colname = 'h2o'
        
        if p_delta == 0:
            self.P_table = init_table.copy()
            self.P_table[water_p_colname] = init_table['hs_p']

        else:
            self.P_table = _integrate_pressure(well_header=well_header,
                                            pt_df_in = init_table, 
                                            get_rho = rho_h2o_getter, 
                                            reference_depth = self.z_resrv, 
                                            reference_pressure = p0, 
                                            direction = 'up', 
                                            colname_p = water_p_colname)
            
            self.P_table = _integrate_pressure(well_header=well_header,
                                            pt_df_in = self.P_table, 
                                            get_rho = rho_h2o_getter, 
                                            reference_depth = self.z_resrv, 
                                            reference_pressure = p0, 
                                            direction = 'down', 
                                            colname_p = water_p_colname)


        
        #CO2
        p0 = np.interp(z_co2_datum, self.P_table['depth_msl'], self.P_table[water_p_colname])
        self.p_CO2_datum = p0
        co2_p_colname  = 'co2'

        self.P_table = _integrate_pressure(well_header=well_header,
                                           pt_df_in = self.P_table, 
                                           get_rho = rho_co2_getter, 
                                           reference_depth = z_co2_datum, 
                                           reference_pressure = p0, 
                                           direction = 'up', 
                                           colname_p = co2_p_colname)

        #We need the density for water given the CO2-pressures, too
        self.P_table = _get_rho_in_pressure_column(self.P_table,
                                                   co2_p_colname, f"h2o_rho_in_co2_column", rho_h2o_getter)
        
        #Compute MSAD
        self.z_MSAD, self.p_MSAD = self.compute_MSAD()

        #Compute delta P
        self.p_delta = self.compute_delta_p(init_table)
        
    def integrate_downwards(self,
                            init_table: pd.DataFrame,
                            well_header: dict,
                            rho_co2_getter: RectBivariateSpline,
                            rho_h2o_getter: RectBivariateSpline,
                            z_co2_datum: float):
        """
        Method to compute maximum pressure tables starting from Shmin at given depth.
        CO2 pressure is integrated downwards from input depth (z_MSAD). 
        H2O pressure is integrated upwards and downwards from point where CO2 pressure meets CO2 datum.
        """
        
        print(f'Pressure scenario {self.name}: Compute maximum pressurization needed to reach Shmin at {self.z_MSAD} mTVDMSL')
        
        #CO2
        co2_p_colname = 'co2'
        self.p_MSAD = np.interp(float(self.z_MSAD), init_table['depth_msl'], init_table['Shmin'])

        self.P_table = _integrate_pressure(well_header=well_header,
                                           pt_df_in = init_table, 
                                           get_rho = rho_co2_getter, 
                                           reference_depth = float(self.z_MSAD), 
                                           reference_pressure = self.p_MSAD, 
                                           direction = 'down', 
                                           colname_p = co2_p_colname)
        
        #Clean up pressure values below Gas_Water_contact
        self.P_table.loc[self.P_table['depth_msl'] > z_co2_datum, co2_p_colname] = np.nan



        #Water
        p0 = np.interp(float(z_co2_datum), self.P_table['depth_msl'], self.P_table[co2_p_colname])
        self.p_CO2_datum = p0
        water_p_colname = 'h2o'

        self.P_table = _integrate_pressure(well_header=well_header,
                                           pt_df_in = self.P_table, 
                                           get_rho = rho_h2o_getter, 
                                           reference_depth = z_co2_datum, 
                                           reference_pressure = p0, 
                                           direction = 'up', 
                                           colname_p = water_p_colname)
        
        self.P_table = _integrate_pressure(well_header=well_header,
                                           pt_df_in = self.P_table, 
                                           get_rho = rho_h2o_getter, 
                                           reference_depth = z_co2_datum, 
                                           reference_pressure = p0, 
                                           direction = 'down', 
                                           colname_p = water_p_colname)

        #We need the density for water given the CO2-pressures, too
        self.P_table = _get_rho_in_pressure_column(self.P_table,
                                                   co2_p_colname, f"h2o_rho_in_co2_column", rho_h2o_getter)

        #Compute reservoir pressure at CO2_datum
        hs_p = self.P_table['h2o'].values

        self.z_resrv = z_co2_datum
        self.p_resrv = self.p_CO2_datum

        #Compute delta P
        self.p_delta = self.compute_delta_p(init_table)

    def compute_MSAD(self):
        """
        Method to compute MSAD.
        """
        depth  = self.P_table['depth_msl'].values
        shmin   = self.P_table['Shmin'].values
        co2_p  = self.P_table['co2'].values
        return compute_intersection(x = depth, y1 = shmin, y2 = co2_p)
            

    def compute_delta_p(self, init_table: pd.DataFrame):
        """
        Method to compute delta P.
        """
        hs_p = np.interp(float(self.z_resrv), init_table['depth_msl'], init_table['hs_p'])
            
        return self.p_resrv - hs_p