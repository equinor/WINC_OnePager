"""
ScenarioParameters: Handles parameter resolution for pressure scenarios.

This class encapsulates the complex logic of resolving scenario parameters
from various combinations of inputs (reservoir pressure, fluid contact, delta pressure, etc.).
"""

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd


@dataclass
class ScenarioParameters:
    """
    Resolves and validates parameters for a pressure scenario.

    Attributes:
        from_resrvr (bool): Whether scenario is computed from reservoir conditions.
        z_fluid_contact (float, optional): Depth of fluid contact.
        p_fluid_contact (float, optional): Pressure at fluid contact.
        z_resrv (float, optional): Reservoir depth.
        p_resrv (float, optional): Reservoir pressure.
        p_delta (float, optional): Pressure difference from hydrostatic.
        z_MSAD (float, optional): Maximum Sustainable Annulus Depth.
        p_MSAD (float, optional): Pressure at MSAD.
        init_curves (pd.DataFrame, optional): Initial curves for interpolation.
    """

    from_resrvr: bool
    z_fluid_contact: Optional[float] = None
    p_fluid_contact: Optional[float] = None
    z_resrv: Optional[float] = None
    p_resrv: Optional[float] = None
    p_delta: Optional[float] = None
    z_MSAD: Optional[float] = None
    p_MSAD: Optional[float] = None
    init_curves: Optional[pd.DataFrame] = None

    def validate_and_resolve(self) -> tuple:
        """
        Validate inputs and resolve parameters to get reference depth and pressure.

        Returns:
            tuple: (reference_depth, reference_pressure)

        Raises:
            ValueError: If insufficient or conflicting parameters are provided.
        """
        if self.from_resrvr:
            return self._resolve_from_reservoir()
        else:
            return self._resolve_from_MSAD()

    def _resolve_from_reservoir(self) -> tuple:
        """
        Resolve parameters when computing from reservoir conditions.

        Returns:
            tuple: (reference_depth, reference_pressure)
        """
        # Create parameter series for easy checking
        params = pd.Series(
            {
                "z_fluid_contact": self.z_fluid_contact,
                "p_fluid_contact": self.p_fluid_contact,
                "p_delta": self.p_delta,
                "p_resrv": self.p_resrv,
                "z_resrv": self.z_resrv,
            },
            dtype=float,
        )

        if params.isna().all():
            raise ValueError("At least one parameter (z_fluid_contact or z_resrv) must be provided.")

        params_with_value = params.dropna()

        # Case 1: Only depth provided (z_fluid_contact or z_resrv) -> hydrostatic
        if params[["p_fluid_contact", "p_resrv", "p_delta"]].isna().all():
            return self._handle_hydrostatic_case(params_with_value)

        # Case 2: z_fluid_contact and p_fluid_contact provided
        if "z_fluid_contact" in params_with_value.index and "p_fluid_contact" in params_with_value.index:
            return self._handle_fluid_contact_case(params_with_value)

        # Case 3: p_delta and z_fluid_contact provided
        if "p_delta" in params_with_value.index and "z_fluid_contact" in params_with_value.index:
            return self._handle_delta_at_contact(params_with_value)

        # Case 4: p_resrv and z_resrv provided
        if "p_resrv" in params_with_value.index and "z_resrv" in params_with_value.index:
            return self.z_resrv, self.p_resrv

        # Case 5: p_resrv and z_fluid_contact provided
        if "p_resrv" in params_with_value.index and "z_fluid_contact" in params_with_value.index:
            self.p_fluid_contact = self.p_resrv
            return self.z_fluid_contact, self.p_fluid_contact

        # Case 6: p_fluid_contact and z_resrv provided
        if "p_fluid_contact" in params_with_value.index and "z_resrv" in params_with_value.index:
            self.z_fluid_contact = self.z_resrv
            return self.z_fluid_contact, self.p_fluid_contact

        # Case 7: p_delta and z_resrv provided
        if "p_delta" in params_with_value.index and "z_resrv" in params_with_value.index:
            self.z_fluid_contact = self.z_resrv
            hydrostatic_at_contact = np.interp(self.z_fluid_contact, self.init_curves["depth"], self.init_curves["hydrostatic_pressure"])
            self.p_fluid_contact = self.p_delta + hydrostatic_at_contact
            return self.z_fluid_contact, self.p_fluid_contact

        raise ValueError("Insufficient or conflicting parameters provided for 'from_resrvr' scenario.")

    def _handle_hydrostatic_case(self, params_with_value: pd.Series) -> tuple:
        """
        Handle case where only depth is provided (hydrostatic assumption).

        Args:
            params_with_value: Series of non-null parameters.

        Returns:
            tuple: (reference_depth, reference_pressure)
        """
        if "z_resrv" in params_with_value.index and self.z_fluid_contact is None:
            self.z_fluid_contact = self.z_resrv
        elif "z_fluid_contact" in params_with_value.index and self.z_resrv is None:
            self.z_resrv = self.z_fluid_contact

        # Get hydrostatic pressure at fluid contact
        self.p_fluid_contact = np.interp(self.z_fluid_contact, self.init_curves["depth"], self.init_curves["hydrostatic_pressure"])
        self.p_delta = 0.0

        return self.z_fluid_contact, self.p_fluid_contact

    def _handle_fluid_contact_case(self, params_with_value: pd.Series) -> tuple:
        """
        Handle case where both z_fluid_contact and p_fluid_contact are provided.

        Args:
            params_with_value: Series of non-null parameters.

        Returns:
            tuple: (reference_depth, reference_pressure)
        """
        if "p_delta" not in params_with_value.index:
            # Compute p_delta
            hydrostatic_at_contact = np.interp(self.z_fluid_contact, self.init_curves["depth"], self.init_curves["hydrostatic_pressure"])
            self.p_delta = self.p_fluid_contact - hydrostatic_at_contact
        else:
            # p_delta also provided - add it to p_fluid_contact
            hydrostatic_at_contact = np.interp(self.z_fluid_contact, self.init_curves["depth"], self.init_curves["hydrostatic_pressure"])
            self.p_fluid_contact = self.p_delta + hydrostatic_at_contact
            self.p_delta = self.p_fluid_contact - hydrostatic_at_contact

        return self.z_fluid_contact, self.p_fluid_contact

    def _handle_delta_at_contact(self, params_with_value: pd.Series) -> tuple:
        """
        Handle case where p_delta and z_fluid_contact are provided.

        Args:
            params_with_value: Series of non-null parameters.

        Returns:
            tuple: (reference_depth, reference_pressure)
        """
        hydrostatic_at_contact = np.interp(self.z_fluid_contact, self.init_curves["depth"], self.init_curves["hydrostatic_pressure"])
        self.p_fluid_contact = self.p_delta + hydrostatic_at_contact

        return self.z_fluid_contact, self.p_fluid_contact

    def _resolve_from_MSAD(self) -> tuple:
        """
        Resolve parameters when computing from MSAD (Maximum Sustainable Annulus Depth).

        Returns:
            tuple: (reference_depth, reference_pressure)
        """
        if self.z_MSAD is None:
            raise ValueError("If 'from_resrvr' is False, you must provide 'z_MSAD'.")

        self.p_MSAD = np.interp(self.z_MSAD, self.init_curves["depth"], self.init_curves["min_horizontal_stress"])

        return self.z_MSAD, self.p_MSAD

    def update_from_profile(self, fluid_pressure_profile: np.ndarray):
        """
        Update derived parameters after fluid pressure profile is computed.

        Args:
            fluid_pressure_profile: Computed fluid pressure curve.
        """
        if not self.from_resrvr:
            return

        depth = self.init_curves["depth"]
        hydrostatic = self.init_curves["hydrostatic_pressure"]

        # Fill in missing parameters based on computed profile
        if self.z_fluid_contact is not None and self.p_fluid_contact is not None:
            if self.p_delta is None:
                self.p_delta = self.p_fluid_contact - np.interp(self.z_fluid_contact, depth, hydrostatic)

            if self.z_resrv is None:
                self.p_resrv = self.p_fluid_contact
                self.z_resrv = self.z_fluid_contact
            elif self.z_resrv < self.z_fluid_contact:
                self.p_resrv = np.interp(self.z_resrv, depth, fluid_pressure_profile)
            else:
                self.p_resrv = self.p_fluid_contact
                self.z_resrv = self.z_fluid_contact

        elif self.p_resrv is not None and self.z_resrv is not None:
            if self.z_fluid_contact is None:
                self.z_fluid_contact = self.z_resrv
                self.p_fluid_contact = self.p_resrv
                self.p_delta = self.p_fluid_contact - np.interp(self.z_fluid_contact, depth, hydrostatic)

            elif self.z_fluid_contact > self.z_resrv:
                self.p_fluid_contact = np.interp(self.z_fluid_contact, depth, fluid_pressure_profile)
                self.p_delta = self.p_fluid_contact - np.interp(self.z_fluid_contact, depth, hydrostatic)
            else:
                self.z_fluid_contact = self.z_resrv
                self.p_fluid_contact = self.p_resrv
                self.p_delta = self.p_fluid_contact - np.interp(self.z_fluid_contact, depth, hydrostatic)

    def update_from_MSAD_intersection(self, z_MSAD: float, p_MSAD: float, z_fluid_contact: float):
        """
        Update parameters when MSAD is found to be deeper than initial fluid contact.

        Args:
            z_MSAD: Computed MSAD depth.
            p_MSAD: Computed MSAD pressure.
            z_fluid_contact: Original fluid contact depth.
        """
        if z_MSAD > z_fluid_contact:
            self.z_MSAD = z_MSAD
            self.p_MSAD = p_MSAD
            self.z_fluid_contact = z_MSAD
            self.p_fluid_contact = p_MSAD
            self.z_resrv = z_MSAD
            self.p_resrv = p_MSAD

            hydrostatic_at_contact = np.interp(self.z_fluid_contact, self.init_curves["depth"], self.init_curves["hydrostatic_pressure"])
            self.p_delta = self.p_fluid_contact - hydrostatic_at_contact
