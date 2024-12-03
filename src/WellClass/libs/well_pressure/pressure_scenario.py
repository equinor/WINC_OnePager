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
    from_resrvr : bool = None
    z_MSAD : float = np.nan
    p_MSAD : float = np.nan
    z_resrv: float = np.nan
    p_resrv: float = np.nan 
    z_fluid_contact: float = np.nan
    p_fluid_contact: float = np.nan
    p_delta: float = np.nan
    use_fluid_sg : bool = False #Use fluid specific gravity from manual input
    fluid_sg  :  float = None  #Manual value of fluid specific gravity    


    def compute_pressure_profile(self,
                                 init_table: pd.DataFrame,
                                 well_header: dict,
                                 rho_co2_getter: RectBivariateSpline,
                                 rho_h2o_getter: RectBivariateSpline
                                 ):
        """
        Depending on input defines which method to run.
        """
        if not isinstance(self.from_resrvr, bool):
            # class will compute max pressure from reservoir upwards
            raise ValueError("Invalid input parameters for computing pressure profile")
        elif self.from_resrvr:
            self.integrate_from_resrv(init_table, well_header,rho_co2_getter, rho_h2o_getter)
        else:
            self.integrate_from_MSAD(init_table, well_header,rho_co2_getter, rho_h2o_getter)
            

    def integrate_from_resrv(self,
                          init_table: pd.DataFrame,
                          well_header: dict,
                          rho_co2_getter: RectBivariateSpline,
                          rho_h2o_getter: RectBivariateSpline):
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
        
        print(f'TEST {p0=}')

        if p_delta == 0:
            self.P_table = init_table.copy()
            self.P_table[water_p_colname] = init_table['hs_p']

        

        else:
            self.P_table = _integrate_pressure(well_header=well_header,
                                            pt_df_in = init_table, 
                                            get_rho = rho_h2o_getter, 
                                            reference_depth = self.z_fluid_contact, 
                                            reference_pressure = p0, 
                                            colname_p = water_p_colname)
            
        
        #CO2
        # print(f'{self.p_resrv=}')
        # p0 = np.interp(self.z_fluid_contact, self.P_table['depth_msl'], self.P_table[water_p_colname])
        # print(f'CO2 datum pressure: {p0:.2f} bar at {self.z_fluid_contact:.2f} mTVDMSL')
        # print(f'CO2 reference pressure: {self.p_resrv:.2f} bar at {self.z_resrv:.2f} mTVDMSL')
        self.p_fluid_contact = p0
        co2_p_colname  = 'co2'

        self.P_table = _integrate_pressure(well_header=well_header,
                                           pt_df_in = self.P_table, 
                                           get_rho = rho_co2_getter, 
                                           reference_depth = self.z_fluid_contact, 
                                           reference_pressure = p0, 
                                           colname_p = co2_p_colname)

        #We need the density for water given the CO2-pressures, too
        self.P_table = _get_rho_in_pressure_column(self.P_table,
                                                   co2_p_colname, f"h2o_rho_in_co2_column", rho_h2o_getter)
        
        #Compute MSAD
        self.z_MSAD, self.p_MSAD = self.compute_MSAD()

        #Compute delta P
        print(f'{self.p_delta=}')
        self.p_delta = self.compute_delta_p(init_table)

        #Clean up pressure values below Gas_Water_contact
        bottom_limit = max(self.z_resrv+1, self.z_fluid_contact)
        self.P_table.loc[self.P_table['depth_msl'] > bottom_limit, co2_p_colname] = np.nan
        self.P_table.loc[self.P_table['depth_msl'] < self.z_MSAD, co2_p_colname] = np.nan
        
    def integrate_from_MSAD(self,
                            init_table: pd.DataFrame,
                            well_header: dict,
                            rho_co2_getter: RectBivariateSpline,
                            rho_h2o_getter: RectBivariateSpline):
        """
        Method to compute maximum pressure tables starting from Shmin at given depth.
        CO2 pressure is integrated downwards from input depth (z_MSAD). 
        H2O pressure is integrated upwards and downwards from point where CO2 pressure meets CO2 datum.
        """
        
        print(f'Pressure scenario {self.name}: Compute maximum pressurization needed to reach Shmin at {self.z_MSAD} mTVDMSL')
        
        #integrate CO2 pressure downwards from MSAD
        co2_p_colname = 'co2'

        # retrieve Shmin (p_MSAD) at MSAD
        self.p_MSAD = np.interp(float(self.z_MSAD), init_table['depth_msl'], init_table['Shmin'])
        print(f'Shmin at MSAD: {self.p_MSAD:.2f} bar at {self.z_MSAD:.2f} mTVDMSL')

        self.P_table = _integrate_pressure(well_header=well_header,
                                           pt_df_in = init_table, 
                                           get_rho = rho_co2_getter, 
                                           reference_depth = float(self.z_MSAD), 
                                           reference_pressure = self.p_MSAD, 
                                           colname_p = co2_p_colname)
        




        #Water
        # retrieve pressure at gas-water contact (p_CO2 == p_water)
        p0 = np.interp(float(self.z_fluid_contact), self.P_table['depth_msl'], self.P_table[co2_p_colname])
        self.p_fluid_contact = p0
        print(f'CO2 datum pressure: {p0:.2f} bar at {self.z_fluid_contact:.2f} mTVDMSL')
        print(f'CO2 datum pressure: {self.p_fluid_contact:.2f} bar at {self.z_fluid_contact:.2f} mTVDMSL')
        water_p_colname = 'h2o'

        self.P_table = _integrate_pressure(well_header=well_header,
                                           pt_df_in = self.P_table, 
                                           get_rho = rho_h2o_getter, 
                                           reference_depth = self.z_fluid_contact, 
                                           reference_pressure = p0, 
                                           colname_p = water_p_colname)
        


        #We need the density for water given the CO2-pressures, too
        self.P_table = _get_rho_in_pressure_column(self.P_table,
                                                   co2_p_colname, f"h2o_rho_in_co2_column", rho_h2o_getter)

        #Compute reservoir pressure at fluid_contact
        self.z_resrv = self.z_fluid_contact
        self.p_resrv = self.p_fluid_contact

        #Compute delta P
        self.p_delta = self.compute_delta_p(init_table)


        #Clean up pressure values below Gas_Water_contact
        bottom_limit = max(self.z_resrv+1, self.z_fluid_contact)
        self.P_table.loc[self.P_table['depth_msl'] > bottom_limit, co2_p_colname] = np.nan
        self.P_table.loc[self.P_table['depth_msl'] < self.z_MSAD, co2_p_colname] = np.nan



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
        hs_p = np.interp(float(self.z_fluid_contact), init_table['depth_msl'], init_table['hs_p'])
            
        return self.p_resrv - hs_p