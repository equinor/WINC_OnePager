from dataclasses import dataclass, field
import pandas as pd
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy import constants as const


from ..pvt.pvt import _integrate_pressure, get_rho_from_pvt_data
from ..utils.compute_intersection import compute_intersection

@dataclass
class PressureScenario:
    name: str
    from_resrvr: bool
    init_curves: pd.DataFrame 
    brine_interpolator: RectBivariateSpline
    fluid_type: str = None
    fluid_interpolator: RectBivariateSpline = None
    fluid_composition: str = None
    pvt_data: dict = None
    specific_gravity: float = None
    z_MSAD: float = None
    p_MSAD: float = None
    z_MSAD_brine: float = None
    p_MSAD_brine: float = None
    p_resrv: float = None
    z_resrv: float = None
    z_fluid_contact: float = None
    p_fluid_contact: float = None
    p_delta: float = None
    cleanup_curves: bool = True

    
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
        self._validate_parameters()
        
        # Compute fluid pressure profile based on provided parameters
        if self.from_resrvr:
            # Handle scenarios when from_resrvr is True
            fluid_pressure_curve = self._handle_from_reservoir()
            self._compute_MSAD(fluid_pressure_curve)
        else:
            # Handle scenarios when from_resrvr is False
            fluid_pressure_curve = self._handle_from_MSAD()

        # Update init_curves table with the computed fluid pressure curve
        self.init_curves['fluid_pressure'] = fluid_pressure_curve
        
        # Compute water pressure profile
        self._compute_brine_pressure() 


        if self.cleanup_curves:
            self._adjust_pressure_curves()

    def _validate_parameters(self):
        """Validate required parameters for pressure computation."""
        if self.name is None:
            raise ValueError("The 'name' parameter is required.")
        if self.from_resrvr is None:
            raise ValueError("The 'from_resrvr' parameter is required.")
        # if self.from_resrvr and self.z_fluid_contact and self.p_delta is None:
        #     self.p_delta = 0


    def _compute_brine_pressure(self):
        """Compute the brine pressure profile based on fluid pressure."""
        print(self.p_delta)
        if np.isclose(self.p_delta, 0, atol=1e-2):
            # If delta_p is zero, use the hydrostatic pressure curve for the water pressure profile
            self.init_curves['brine_pressure'] = self.init_curves['hydrostatic_pressure']
        else:
            # If delta_p is not zero, integrate using z_fluid_contact and p_fluid_contact
            self.init_curves['brine_pressure'] = self._integrate_brine_pressure_curve( reference_depth=self.z_fluid_contact,
                                                                                     reference_pressure=self.p_fluid_contact)

            if self.p_delta > 0:
                self.z_MSAD_brine, self.p_MSAD_brine = compute_intersection(self.init_curves['depth'].values, self.init_curves['brine_pressure'].values, self.init_curves['min_horizontal_stress'].values)


    def _handle_from_reservoir(self) -> np.ndarray:
        """ 
        Method to handle scenarios when the fluid pressure profile is computed from the reservoir.
        """

        ip_params = zip(
            ['z_fluid_contact', 'p_fluid_contact', 'p_delta', 'p_resrv', 'z_resrv'],
            np.array([self.z_fluid_contact, self.p_fluid_contact, self.p_delta, self.p_resrv, self.z_resrv], dtype=float)
        )

        ip_params = pd.Series(dict(ip_params))


        # Check if only z_fluid_contact is provided (hydrostatic default case)
        # if all(param is None for param in (self.p_delta, self.p_resrv, self.z_resrv, self.p_fluid_contact, self.z_fluid_contact)):
        if ip_params.isna().all():
            raise ValueError("At least one parameter (z_fluid_contact or z_resrv) must be provided.")
        

        ip_params_with_value = ip_params.dropna()


        # Check if only z_fluid_contact or z_resrv are provided, make a hydrostatic case
        if ip_params[['p_fluid_contact', 'p_resrv', 'p_delta']].isna().all() and ('z_fluid_contact' in ip_params_with_value.index or 'z_resrv' in ip_params_with_value.index):
            if self.z_resrv is None:
                # Only z_fluid_contact is provided, assume hydrostatic pressure
                self.z_resrv = self.z_fluid_contact

            if self.z_fluid_contact is None:
                self.z_fluid_contact = self.z_resrv
            
            
            self.p_fluid_contact = np.interp(self.z_fluid_contact, self.init_curves['depth'], self.init_curves['hydrostatic_pressure'])

            ref_z = self.z_fluid_contact
            ref_p = self.p_fluid_contact

        
        elif 'z_fluid_contact' in ip_params_with_value.index and 'p_fluid_contact' in ip_params_with_value.index:
            # Both z_fluid_contact and p_fluid_contact are provided
            self.p_delta = self.p_fluid_contact - np.interp(self.z_fluid_contact, self.init_curves['depth'], self.init_curves['hydrostatic_pressure'])

            ref_z = self.z_fluid_contact
            ref_p = self.p_fluid_contact

            
        
        elif 'p_delta' in ip_params_with_value.index and 'z_fluid_contact' in ip_params_with_value.index:
            # Both p_delta and z_fluid_contact are provided
            self.p_fluid_contact = self.p_delta + np.interp(self.z_fluid_contact, self.init_curves['depth'], self.init_curves['hydrostatic_pressure'])

            ref_z = self.z_fluid_contact
            ref_p = self.p_fluid_contact

        elif 'p_resrv' in ip_params_with_value.index and 'z_resrv' in ip_params_with_value.index:
            ref_z = self.z_resrv
            ref_p = self.p_resrv

        elif 'p_resrv' in ip_params_with_value.index and 'z_fluid_contact' in ip_params_with_value.index:
            self.p_fluid_contact = self.p_resrv
            ref_z = self.z_fluid_contact
            ref_p = self.p_fluid_contact

        elif 'p_fluid_contact' in ip_params_with_value.index and 'z_resrv' in ip_params_with_value.index:
            self.z_fluid_contact = self.z_resrv
            ref_z = self.z_fluid_contact
            ref_p = self.p_fluid_contact

        elif 'p_delta' in ip_params_with_value.index and 'z_resrv' in ip_params_with_value.index:
            self.z_fluid_contact = self.z_resrv
            self.p_fluid_contact = self.p_delta + np.interp(self.z_fluid_contact, self.init_curves['depth'], self.init_curves['hydrostatic_pressure'])
            ref_z = self.z_fluid_contact
            ref_p = self.p_fluid_contact

        else:
            ValueError("Insufficient parameters provided for 'from_resrvr' scenario.")            
        

        


        # print('\nafter check\n', ip_params)

        fluid_pressure_profile = self._compute_fluid_pressure_curve(
            reference_depth=ref_z,
            reference_pressure=ref_p,
            fluid_key=self.fluid_type
        )


        # Fill missing parameters based on combination of provided parameters
        if self.z_fluid_contact and self.p_fluid_contact:
            if self.p_delta is None:
            # Compute p_delta if z_fluid_contact and p_fluid_contact are provided
                self.p_delta = self.p_fluid_contact - np.interp(self.z_fluid_contact, self.init_curves['depth'], self.init_curves['hydrostatic_pressure'])
            
            if self.z_resrv is None:
                self.p_resrv = self.p_fluid_contact
                self.z_resrv = self.z_fluid_contact

            elif self.z_resrv < self.z_fluid_contact:
                # If z_resrv is less than z_fluid_contact, set p_resrv to the fluid pressure at z_fluid_contact
                self.p_resrv = np.interp(self.z_resrv, self.init_curves['depth'], fluid_pressure_profile)
            
            else:
                # If z_resrv is greater than or equal to z_fluid_contact, set p_resrv to p_fluid_contact
                self.p_resrv = self.p_fluid_contact
                self.z_resrv = self.z_fluid_contact

        elif self.p_resrv and self.z_resrv:
            if self.z_fluid_contact is None:
                # If z_fluid_contact is not provided, set it to z_resrv
                self.z_fluid_contact = self.z_resrv
                self.p_fluid_contact = self.p_resrv
                self.p_delta = self.p_fluid_contact - np.interp(self.z_fluid_contact, self.init_curves['depth'], self.init_curves['hydrostatic_pressure'])

            elif self.z_fluid_contact > self.z_resrv:
                # If z_fluid_contact is greater than z_resrv, set p_fluid_contact to p_resrv
                self.p_fluid_contact = np.interp(self.z_fluid_contact, self.init_curves['depth'], fluid_pressure_profile)
                self.p_delta = self.p_fluid_contact - np.interp(self.z_fluid_contact, self.init_curves['depth'], self.init_curves['hydrostatic_pressure'])

            else:
                # If z_fluid_contact is less than or equal to z_resrv, set p_fluid_contact to p_resrv
                self.z_fluid_contact = self.z_resrv
                self.p_fluid_contact = self.p_resrv
                self.p_delta = self.p_fluid_contact - np.interp(self.z_fluid_contact, self.init_curves['depth'], self.init_curves['hydrostatic_pressure'])
            
        

        #PRINT DEBUGGING INFO
        # op_params = dict(zip(
        #     ['z_fluid_contact', 'p_fluid_contact', 'p_delta', 'p_resrv', 'z_resrv'],
        #     np.array([self.z_fluid_contact, self.p_fluid_contact, self.p_delta, self.p_resrv, self.z_resrv], dtype=float)
        # ))

        
        # for key, ip_param, op_param in zip(ip_params.keys(), ip_params, pd.Series(op_params)):
        #     if np.isclose(ip_param, op_param):
        #         # If the input parameter is close to the output parameter, update the output parameter
        #         print(f'{key}: Input value used.')
        #     elif np.isnan(ip_param):
        #         # If the input parameter is None, set it to the output parameter
        #         print(f'{key}: value computed and updated.')
        #     elif ~np.isnan(ip_param):
        #         # If the input parameter is not None, keep the input value
        #         print(f'{key}: Input value {ip_param} overriden by {op_param}.')
        return fluid_pressure_profile



    def _compute_MSAD(self, fluid_pressure_profile):

        if 'min_horizontal_stress' in self.init_curves.columns:
            depth = self.init_curves['depth'].values
            shmin = self.init_curves['min_horizontal_stress'].values

            self.z_MSAD, self.p_MSAD = compute_intersection(depth, fluid_pressure_profile, shmin)



            # if self.p_delta > 0:
            #     self.z_MSAD_brine, self.p_MSAD_brine = compute_intersection(depth, self.init_curves['brine_pressure'].values, shmin)


            if self.z_MSAD > self.z_fluid_contact:
                # If z_MSAD is greater than or equal to z_fluid_contact, set them to NaN
                self.z_fluid_contact = self.z_MSAD
                self.p_fluid_contact = self.p_MSAD
                self.z_resrv = self.z_MSAD
                self.p_resrv = self.p_MSAD
                self.p_delta = self.p_fluid_contact - np.interp(self.z_fluid_contact, self.init_curves['depth'], self.init_curves['hydrostatic_pressure'])

        else:
            raise KeyError("Shmin column is missing from init_curves DataFrame")

        

    def _handle_from_MSAD(self) -> np.ndarray:

        if self.z_MSAD is None:
            raise ValueError("If 'from_resrvr' is False, you must provide 'z_MSAD'.")
        
        self.p_MSAD = np.interp(self.z_MSAD, self.init_curves['depth'], self.init_curves['min_horizontal_stress'])

        fluid_pressure_curve = self._compute_fluid_pressure_curve( reference_depth=self.z_MSAD,
                                                                                reference_pressure=self.p_MSAD,
                                                                                fluid_key=self.fluid_type)
        
        # Check if z_fluid_contact is provided
        

        if self.z_fluid_contact is not None and  self.z_MSAD > self.z_fluid_contact:
            # z_MSAD to evaluate deeper than z_fluid_contact.
            self.z_fluid_contact = self.z_MSAD
            self.p_fluid_contact = self.p_MSAD
            self.z_resrv = self.z_MSAD
            self.p_resrv = self.p_MSAD
            self.p_delta = self.p_fluid_contact - np.interp(self.z_fluid_contact, self.init_curves['depth'], self.init_curves['hydrostatic_pressure'])
            
        
        
        elif self.z_fluid_contact is not None:
            self.p_fluid_contact = np.interp(self.z_fluid_contact, self.init_curves['depth'], fluid_pressure_curve)
            self.p_delta = self.p_fluid_contact - np.interp(self.z_fluid_contact, self.init_curves['depth'], self.init_curves['hydrostatic_pressure'])



        else:
            # Assumes hydrostatic pressure at reservoir
            self.p_delta = 0
            self.z_fluid_contact, self.p_fluid_contact = compute_intersection(self.init_curves['depth'].values, fluid_pressure_curve, self.init_curves['hydrostatic_pressure'].values)
            
        self.p_resrv = self.p_fluid_contact
        self.z_resrv = self.z_fluid_contact
        
        return fluid_pressure_curve
        
    def _adjust_pressure_curves(self):
        """
        Adjust the fluid_pressure_curve and brine_pressure_curve to be valid within the given ranges.
        """


        for z_val in set([self.z_MSAD, self.z_fluid_contact, self.z_resrv]):

            # Check if the z_val is close to any depth in the init_curves
            eval_array = ~np.isclose(self.init_curves['depth'], z_val, atol=0.001, rtol=0)
            if np.all(eval_array) and ~np.isnan(z_val):

                new_index = z_val
                new_record = {'depth': new_index }

                # insert the new record into the DataFrame
                self.init_curves.loc[new_index] = new_record

        self.init_curves = self.init_curves.sort_index()
        self.init_curves = self.init_curves.interpolate()

        self.init_curves.loc[self.init_curves['fluid_pressure'] > self.init_curves['min_horizontal_stress'], 'fluid_pressure'] = np.nan
        self.init_curves.loc[self.init_curves['brine_pressure'] > self.init_curves['min_horizontal_stress'], 'brine_pressure'] = np.nan
        self.init_curves.loc[self.init_curves['fluid_pressure'] < self.init_curves['brine_pressure'], 'fluid_pressure'] = np.nan




    def _integrate_brine_pressure_curve(self, 
                                      reference_depth: float, 
                                      reference_pressure: float) -> np.ndarray:

        return _integrate_pressure(
                init_curves=self.init_curves,
                reference_depth=reference_depth,
                reference_pressure=reference_pressure,
                pvt_data=self.pvt_data,
                fluid_key='brine',
                interpolator=self.brine_interpolator
            )


    def _compute_fluid_pressure_curve(self, 
                                      reference_depth: float, 
                                      reference_pressure: float,
                                      fluid_key: str = None) -> np.ndarray:
        """Computes the fluid pressure curve based on a reference depth and pressure.

        This method can be used to compute the pressure curve for the specified fluid type
        or for brine if 'brine' is passed as the fluid_key.

        Args:
            reference_depth (float): The depth at which the reference pressure is known.
            reference_pressure (float): The known pressure at the reference depth.
            fluid_key (str, optional): Key indicating the type of fluid or 'brine'. 
                                    Defaults to the scenario's fluid type.

        Returns:
            np.ndarray: The computed fluid pressure curve.
        """
        depth_array = self.init_curves['depth'].values

        if self.specific_gravity is not None:
            # Convert specific gravity to density (assuming specific gravity is relative to water at 4°C, 1000 kg/m³)
            density = self.specific_gravity * 1000  # kg/m³
            # Use the constant density to compute the pressure curve directly
            pressure_curve = reference_pressure + (depth_array - reference_depth) * density * const.g / const.bar
        else:
            # Use PVT data and integration to compute the pressure curve
            interpolator = self.fluid_interpolator
            pressure_curve = _integrate_pressure(
                init_curves=self.init_curves,
                reference_depth=reference_depth,
                reference_pressure=reference_pressure,
                pvt_data=self.pvt_data,
                fluid_key=self.fluid_type,
                interpolator=interpolator
            )
        return pressure_curve