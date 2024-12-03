from dataclasses import dataclass, field
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

from .new_aux_func import _integrate_pressure, get_rho_from_pvt_data
from ..utils.compute_intersection import compute_intersection

@dataclass
class PressureScenario:
    name: str
    fluid_type: str
    pvt_data: dict
    from_resrvr: bool
    init_curves: pd.DataFrame 
    z_MSAD: float = None
    p_MSAD: float = None
    p_resrv: float = None
    z_resrv: float = None
    z_fluid_contact: float = None
    p_fluid_contact: float = None
    p_delta: float = None
    
    def compute_pressure_profile(self):
        print(f'Computing pressure profile for scenario: {self.name}')
        # Create default NaN arrays for fluid and brine pressure
        default_length = len(self.init_curves['depth'])
        default_nan_array = np.full(default_length, np.nan)

        # Begin with assumption that fluid and brine pressure profiles
        # will be set to NaNs if the inputs are insufficient or incorrect
        self.init_curves['fluid_pressure'] = default_nan_array
        self.init_curves['brine_pressure'] = default_nan_array

        # Validate required parameters
        if self.name is None:
            # Return error ir scenario name is not provided
            raise ValueError("The 'name' parameter is required.")
        
        if self.from_resrvr is None:
            # Return error if 'from_resrvr' is not provided
            raise ValueError("The 'from_resrvr' parameter is required.")
        
        # Compute fluid pressure profile based on provided parameters
        if self.from_resrvr:
            # Handle scenarios when from_resrvr is True
            print(f"From Reservoir {self.name=}")
            fluid_pressure_curve = self._handle_from_reservoir()
        else:
            # Handle scenarios when from_resrvr is False
            fluid_pressure_curve = self._handle_from_MSAD()

        # Update init_curves table with the computed fluid pressure curve
        self.init_curves['fluid_pressure'] = fluid_pressure_curve
        
        
        # Compute water pressure profile
        if self.p_delta == 0:
            # If delta_p is zero, use the hydrostatic pressure curve for the water pressure profile
            self.init_curves['brine_pressure'] = self.init_curves['hydrostatic_pressure']
        else:
            # If delta_p is not zero, integrate using z_fluid_contact and p_fluid_contact
            self.init_curves['brine_pressure'] = self._compute_fluid_pressure_curve( reference_depth=self.z_fluid_contact,
                                                                                     reference_pressure=self.p_fluid_contact,
                                                                                     fluid_key='brine')            


    def _handle_from_reservoir(self) -> np.ndarray:
        """ 
        Method to handle scenarios when the fluid pressure profile is computed from the reservoir.
        """

        # Check if only z_fluid_contact is provided (hydrostatic default case)
        all_other_none = all(x is None for x in [self.p_delta, self.p_resrv, self.z_resrv, self.p_fluid_contact])
        print(f'{all_other_none=}')
        for var_name in ['p_delta', 'p_resrv', 'z_resrv', 'p_fluid_contact', 'z_fluid_contact']:    
            value = getattr(self, var_name)
            print(f'{var_name}={value}')
        
        if self.z_fluid_contact is not None and all_other_none:
            # Assume delta_p is zero
            self.p_delta = 0
            # Compute hydrostatic pressure at z_fluid_contact
            self.p_fluid_contact = np.interp(self.z_fluid_contact, self.init_curves['depth'], self.init_curves['hydrostatic_pressure'])
            # Set reservoir pressure and depth to fluid contact values
            self.p_resrv = self.p_fluid_contact
            self.z_resrv = self.z_fluid_contact
            # Compute fluid pressure profile from z_fluid_contact
            fluid_pressure_profile = self._compute_fluid_pressure_curve(
                reference_depth=self.z_fluid_contact,
                reference_pressure=self.p_fluid_contact
            )
        
        # Check if z_fluid_contact and p_delta are provided
        elif self.p_delta is not None and self.z_fluid_contact is not None:
            print('TESTING')
            # Compute fluid pressure profile starting from z_fluid_contact
            self.p_fluid_contact = self.p_delta + np.interp(self.z_fluid_contact, self.init_curves['depth'], self.init_curves['hydrostatic_pressure'])
            fluid_pressure_profile = self._compute_fluid_pressure_curve( reference_depth=self.z_fluid_contact,
                                                                            reference_pressure=self.p_fluid_contact)
        
            if self.z_resrv is None:
                self.z_resrv = self.z_fluid_contact
                self.p_resrv = self.p_fluid_contact
            else:
                if self.z_resrv < self.z_fluid_contact:
                    self.p_resrv = np.interp(self.z_resrv, self.init_curves['depth'], fluid_pressure_profile)
                # else: ignore (z_resrv >= z_fluid_contact)

        # Check if p_resrv and z_resrv are provided
        elif self.p_resrv is not None and self.z_resrv is not None:
            # Compute fluid pressure profile starting from z_resrv
            fluid_pressure_profile = self._compute_fluid_pressure_curve( reference_depth=self.z_resrv,
                                                                    reference_pressure=self.p_resrv)

            if self.z_fluid_contact is None:
                self.p_fluid_contact = self.p_resrv
                self.z_fluid_contact = self.z_resrv
                self.p_delta = self.p_fluid_contact - np.interp(self.z_fluid_contact, self.init_curves['depth'], self.init_curves['hydrostatic_pressure'])
            else:
                if self.z_resrv < self.z_fluid_contact:
                    self.p_fluid_contact = np.interp(self.z_fluid_contact, self.init_curves['depth'], fluid_pressure_profile)
                else:
                    self.p_fluid_contact = self.p_resrv
                    self.z_fluid_contact = self.z_resrv
                    self.p_delta = self.p_fluid_contact - np.interp(self.z_fluid_contact, self.init_curves['depth'], self.init_curves['hydrostatic_pressure'])
                
        else:
            raise ValueError("Insufficient parameters provided for 'from_resrvr' scenario.")

        if 'min_horizontal_stress' in self.init_curves.columns:
            depth = self.init_curves['depth'].values
            shmin = self.init_curves['min_horizontal_stress'].values

            self.z_MSAD, self.p_MSAD = compute_intersection(depth, fluid_pressure_profile, shmin)
        else:
            raise KeyError("Shmin column is missing from init_curves DataFrame")

        return fluid_pressure_profile
        

    def _handle_from_MSAD(self) -> np.ndarray:

        if self.z_MSAD is None:
            raise ValueError("If 'from_resrvr' is False, you must provide 'z_MSAD'.")
        
        self.p_MSAD = np.interp(self.z_MSAD, self.init_curves['depth'], self.init_curves['min_horizontal_stress'])

        fluid_pressure_curve = self._compute_fluid_pressure_curve( reference_depth=self.z_MSAD,
                                                                                reference_pressure=self.p_MSAD)
        
        # Check if z_fluid_contact is provided
        if self.z_fluid_contact is not None:
            self.p_fluid_contact = np.interp(self.z_fluid_contact, self.init_curves['depth'], fluid_pressure_curve)
            self.p_delta = self.p_fluid_contact - np.interp(self.z_fluid_contact, self.init_curves['depth'], self.init_curves['hydrostatic_pressure'])



        else:
            # Assumes hydrostatic pressure at reservoir
            self.p_delta = 0
            self.z_fluid_contact, self.p_fluid_contact = compute_intersection(self.init_curves['depth'].values, fluid_pressure_curve, self.init_curves['hydrostatic_pressure'].values)
            
        self.p_resrv = self.p_fluid_contact
        self.z_resrv = self.z_fluid_contact
        
        return fluid_pressure_curve
        

    def _compute_fluid_pressure_curve(self, 
                                      reference_depth: float, 
                                      reference_pressure: float,
                                      fluid_key: str = None) -> np.ndarray:
        # Wrapper method to call _integrate_pressure with the necessary arguments
        
        
        if fluid_key is None:
            fluid_key = self.fluid_type

       
        return _integrate_pressure(
            init_curves=self.init_curves,
            reference_depth=reference_depth,
            reference_pressure=reference_pressure,
            pvt_data=self.pvt_data,
            fluid_key=fluid_key,
            get_rho_func=get_rho_from_pvt_data
        )
