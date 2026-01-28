from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import constants as const
from scipy.interpolate import RectBivariateSpline

from ..pvt.pvt import _integrate_pressure, get_rho_from_pvt_data
from ..utils.compute_intersection import compute_intersection
from .ScenarioParameters import ScenarioParameters


@dataclass
class PressureScenario:
    """
    Represents a single pressure scenario with fluid and brine pressure profiles.

    This class computes pressure profiles for a specific scenario based on
    reservoir conditions or MSAD (Maximum Sustainable Annulus Depth).
    """

    name: str
    from_resrvr: bool
    init_curves: pd.DataFrame
    brine_interpolator: RectBivariateSpline
    fluid_type: str = None
    fluid_interpolator: RectBivariateSpline = None
    fluid_composition: str = None
    pvt_data: dict = None
    specific_gravity: float = None
    rho_brine: float = None
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
    warnings: list = None

    def compute_pressure_profile(self):
        """Compute fluid and brine pressure profiles for this scenario."""
        print(f"Computing pressure profile for scenario: {self.name}")

        # Create default NaN arrays for fluid and brine pressure
        default_length = len(self.init_curves["depth"])
        default_nan_array = np.full(default_length, np.nan)

        # Begin with assumption that fluid and brine pressure profiles
        # will be set to NaNs if the inputs are insufficient or incorrect
        self.init_curves["fluid_pressure"] = default_nan_array
        self.init_curves["brine_pressure"] = default_nan_array

        # Use ScenarioParameters to resolve reference depth and pressure
        params = ScenarioParameters(
            from_resrvr=self.from_resrvr,
            z_fluid_contact=self.z_fluid_contact,
            p_fluid_contact=self.p_fluid_contact,
            z_resrv=self.z_resrv,
            p_resrv=self.p_resrv,
            p_delta=self.p_delta,
            z_MSAD=self.z_MSAD,
            p_MSAD=self.p_MSAD,
            init_curves=self.init_curves,
        )

        ref_depth, ref_pressure = params.validate_and_resolve()

        # Compute fluid pressure profile
        fluid_pressure_curve = self._compute_fluid_pressure_curve(
            reference_depth=ref_depth, reference_pressure=ref_pressure, fluid_key=self.fluid_type
        )

        # Update parameters from computed profile
        params.update_from_profile(fluid_pressure_curve)

        # Update instance attributes from resolved parameters
        self._update_from_params(params)

        # Compute MSAD if coming from reservoir
        if self.from_resrvr:
            self._compute_MSAD(fluid_pressure_curve)
        else:
            # For MSAD scenarios, ensure p_delta is calculated
            if self.z_fluid_contact is not None:
                if self.p_delta is None:
                    self.p_fluid_contact = np.interp(self.z_fluid_contact, self.init_curves["depth"], fluid_pressure_curve)
                    self.p_delta = self.p_fluid_contact - np.interp(
                        self.z_fluid_contact, self.init_curves["depth"], self.init_curves["hydrostatic_pressure"]
                    )
            else:
                # If no fluid contact specified, set p_delta to 0
                if self.p_delta is None:
                    self.p_delta = 0.0

        # Update init_curves table with the computed fluid pressure curve
        self.init_curves["fluid_pressure"] = fluid_pressure_curve

        # Compute water pressure profile
        self._compute_brine_pressure()

        if self.cleanup_curves:
            self._adjust_pressure_curves()

    def _update_from_params(self, params: ScenarioParameters):
        """Update instance attributes from resolved parameters."""
        self.z_fluid_contact = params.z_fluid_contact
        self.p_fluid_contact = params.p_fluid_contact
        self.z_resrv = params.z_resrv
        self.p_resrv = params.p_resrv
        self.p_delta = params.p_delta
        self.z_MSAD = params.z_MSAD
        self.p_MSAD = params.p_MSAD

    def _compute_brine_pressure(self):
        """Compute the brine pressure profile based on fluid pressure."""
        # Handle None p_delta (shouldn't happen but defensive programming)
        if self.p_delta is None:
            self.p_delta = 0.0

        if np.isclose(self.p_delta, 0, atol=1e-2):
            # If delta_p is zero, use the hydrostatic pressure curve for the water pressure profile
            self.init_curves["brine_pressure"] = self.init_curves["hydrostatic_pressure"]

        elif self.rho_brine is not None:
            self.init_curves["brine_pressure"] = self.init_curves["hydrostatic_pressure"] + self.p_delta

            if self.p_delta > 0:
                self.z_MSAD_brine, self.p_MSAD_brine = compute_intersection(
                    self.init_curves["depth"].values, self.init_curves["brine_pressure"].values, self.init_curves["min_horizontal_stress"].values
                )

        else:
            # If delta_p is not zero, integrate using z_fluid_contact and p_fluid_contact
            brine_pressure, _ = self._integrate_brine_pressure_curve(reference_depth=self.z_fluid_contact, reference_pressure=self.p_fluid_contact)
            self.init_curves["brine_pressure"] = brine_pressure

            if self.p_delta > 0:
                self.z_MSAD_brine, self.p_MSAD_brine = compute_intersection(
                    self.init_curves["depth"].values, self.init_curves["brine_pressure"].values, self.init_curves["min_horizontal_stress"].values
                )

    def _compute_MSAD(self, fluid_pressure_profile: np.ndarray):
        """Compute Maximum Sustainable Annulus Depth (MSAD) from fluid pressure and Shmin intersection."""
        if "min_horizontal_stress" not in self.init_curves.columns:
            raise KeyError("Shmin column is missing from init_curves DataFrame")

        depth = self.init_curves["depth"].values
        shmin = self.init_curves["min_horizontal_stress"].values

        self.z_MSAD, self.p_MSAD = compute_intersection(depth, fluid_pressure_profile, shmin)

        # Update fluid contact if MSAD is deeper or if fluid contact is not set
        if self.z_fluid_contact is None or self.z_MSAD > self.z_fluid_contact:
            # If z_MSAD is greater than z_fluid_contact, update all contact parameters
            self.z_fluid_contact = self.z_MSAD
            self.p_fluid_contact = self.p_MSAD
            self.z_resrv = self.z_MSAD
            self.p_resrv = self.p_MSAD
            self.p_delta = self.p_fluid_contact - np.interp(self.z_fluid_contact, self.init_curves["depth"], self.init_curves["hydrostatic_pressure"])

    def _adjust_pressure_curves(self):
        """
        Adjust the fluid_pressure_curve and brine_pressure_curve to be valid within the given ranges.
        """

        for z_val in set([self.z_MSAD, self.z_fluid_contact, self.z_resrv]):
            # Skip None values
            if z_val is None or np.isnan(z_val):
                continue
            # Check if the z_val is close to any depth in the init_curves
            eval_array = ~np.isclose(self.init_curves["depth"], z_val, atol=0.001, rtol=0)
            if np.all(eval_array):
                new_index = z_val
                new_record = {"depth": new_index}

                # insert the new record into the DataFrame
                self.init_curves.loc[new_index] = new_record

        self.init_curves = self.init_curves.sort_index()
        self.init_curves = self.init_curves.interpolate()

        self.init_curves.loc[self.init_curves["fluid_pressure"] > self.init_curves["min_horizontal_stress"], "fluid_pressure"] = np.nan
        self.init_curves.loc[self.init_curves["brine_pressure"] > self.init_curves["min_horizontal_stress"], "brine_pressure"] = np.nan
        self.init_curves.loc[self.init_curves["fluid_pressure"] < self.init_curves["brine_pressure"], "fluid_pressure"] = np.nan

    def _integrate_brine_pressure_curve(self, reference_depth: float, reference_pressure: float) -> np.ndarray:
        return _integrate_pressure(
            init_curves=self.init_curves,
            reference_depth=reference_depth,
            reference_pressure=reference_pressure,
            pvt_data=self.pvt_data,
            fluid_key="brine",
            interpolator=self.brine_interpolator,
        )

    def _compute_fluid_pressure_curve(self, reference_depth: float, reference_pressure: float, fluid_key: str = None) -> np.ndarray:
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
        depth_array = self.init_curves["depth"].values

        if self.specific_gravity is not None:
            # Convert specific gravity to density (assuming specific gravity is relative to water at 4°C, 1000 kg/m³)
            density = self.specific_gravity * 1000  # kg/m³
            # Use the constant density to compute the pressure curve directly
            pressure_curve = reference_pressure + (depth_array - reference_depth) * density * const.g / const.bar
        else:
            # Use PVT data and integration to compute the pressure curve
            interpolator = self.fluid_interpolator
            pressure_curve, warnings = _integrate_pressure(
                init_curves=self.init_curves,
                reference_depth=reference_depth,
                reference_pressure=reference_pressure,
                pvt_data=self.pvt_data,
                fluid_key=self.fluid_type,
                interpolator=interpolator,
            )

            self.warnings = warnings

        return pressure_curve
