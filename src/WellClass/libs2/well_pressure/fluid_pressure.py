from dataclasses import dataclass, field

import numpy as np
from scipy import constants as const
from scipy.interpolate import interp1d

from .pressure_table import PressureTable


@dataclass
class FluidPressure:
    """
    A single fluid-pressure scenario computed against a shared PressureTable.

    The user provides a reference point (z_fluid_contact + p_fluid_contact,
    or z_fluid_contact + p_delta, or z_resrv + p_resrv) and a fluid density.
    The class computes a constant-density fluid column from that reference
    and derives brine overpressure accordingly.
    """

    # Shared background curves
    table: PressureTable

    # Scenario metadata
    name: str = "default"

    # Reference point parameters (at least one pair must be provided)
    z_fluid_contact: float | None = None
    p_fluid_contact: float | None = None  # pressure at fluid contact (bar)
    p_delta: float | None = None  # overpressure above hydrostatic (bar)
    z_resrv: float | None = None
    p_resrv: float | None = None  # reservoir pressure (bar)

    # Fluid density
    rho_fluid: float = 400.0  # fluid density (kg/m³), e.g. CO2

    # Computed arrays
    fluid_pressure: np.ndarray = field(init=False)
    brine_pressure: np.ndarray = field(init=False)

    # Computed scalars
    z_msad: float | None = field(init=False, default=None)
    p_msad: float | None = field(init=False, default=None)

    def __post_init__(self) -> None:
        self._resolve_reference_point()
        self.fluid_pressure = self._compute_fluid_pressure()
        self.brine_pressure = self._compute_brine_pressure()
        self.z_msad, self.p_msad = self._compute_msad()

    # --- convenience accessors for the shared table arrays ---
    @property
    def depth(self) -> np.ndarray:
        return self.table.depth

    @property
    def hydrostatic_pressure(self) -> np.ndarray:
        return self.table.hydrostatic_pressure

    @property
    def min_horizontal_stress(self) -> np.ndarray:
        return self.table.min_horizontal_stress

    # ------------------------------------------------------------------
    # Reference point resolution
    # ------------------------------------------------------------------
    def _resolve_reference_point(self) -> None:
        """Derive z_fluid_contact, p_fluid_contact and p_delta from whatever the user provided."""
        if self.z_fluid_contact is not None and self.p_fluid_contact is not None:
            hydro_at_contact = float(np.interp(self.z_fluid_contact, self.depth, self.hydrostatic_pressure))
            if self.p_delta is None:
                self.p_delta = self.p_fluid_contact - hydro_at_contact

        elif self.z_fluid_contact is not None and self.p_delta is not None:
            hydro_at_contact = float(np.interp(self.z_fluid_contact, self.depth, self.hydrostatic_pressure))
            self.p_fluid_contact = hydro_at_contact + self.p_delta

        elif self.z_resrv is not None and self.p_resrv is not None:
            self.z_fluid_contact = self.z_resrv
            self.p_fluid_contact = self.p_resrv
            hydro_at_contact = float(np.interp(self.z_fluid_contact, self.depth, self.hydrostatic_pressure))
            self.p_delta = self.p_fluid_contact - hydro_at_contact

        elif self.z_fluid_contact is not None:
            # Hydrostatic default – no overpressure
            self.p_fluid_contact = float(np.interp(self.z_fluid_contact, self.depth, self.hydrostatic_pressure))
            self.p_delta = 0.0

        else:
            raise ValueError(
                "Provide at least (z_fluid_contact) or (z_fluid_contact + p_fluid_contact) or (z_fluid_contact + p_delta) or (z_resrv + p_resrv)."
            )

    # ------------------------------------------------------------------
    # Fluid pressure
    # ------------------------------------------------------------------
    def _compute_fluid_pressure(self) -> np.ndarray:
        """Compute fluid pressure using constant-density column from the reference point."""
        # P(z) = P_ref + rho_fluid * g * (z - z_ref)   [Pa -> bar]
        delta_depth = self.depth - self.z_fluid_contact
        fluid_p = self.p_fluid_contact + (self.rho_fluid * const.g * delta_depth) / 1e5
        return fluid_p

    # ------------------------------------------------------------------
    # Brine pressure
    # ------------------------------------------------------------------
    def _compute_brine_pressure(self) -> np.ndarray:
        """Compute brine (water) pressure as hydrostatic + overpressure."""
        if np.isclose(self.p_delta, 0.0, atol=1e-2):
            return self.hydrostatic_pressure.copy()
        return self.hydrostatic_pressure + self.p_delta

    # ------------------------------------------------------------------
    # MSAD – intersection of fluid pressure with Shmin
    # ------------------------------------------------------------------
    def _compute_msad(self) -> tuple[float | None, float | None]:
        """Find depth where fluid pressure meets minimum horizontal stress."""
        diff = self.fluid_pressure - self.min_horizontal_stress
        sign_changes = np.where(np.diff(np.sign(diff)))[0]

        if len(sign_changes) == 0:
            return None, None

        # Take the shallowest crossing
        idx = sign_changes[0]
        z_msad = float(np.interp(0, [diff[idx], diff[idx + 1]], [self.depth[idx], self.depth[idx + 1]]))
        p_msad = float(np.interp(z_msad, self.depth, self.min_horizontal_stress))
        return z_msad, p_msad

    # ------------------------------------------------------------------
    # Query
    # ------------------------------------------------------------------
    def get_values_at_depth(self, depth_value: float) -> dict[str, float]:
        """Interpolate all curves at a given depth."""
        base = self.table.get_values_at_depth(depth_value)

        def _interp(arr: np.ndarray) -> float:
            f = interp1d(self.depth, arr, bounds_error=False, fill_value="extrapolate")
            return float(f(depth_value))

        base["fluid_pressure"] = _interp(self.fluid_pressure)
        base["brine_pressure"] = _interp(self.brine_pressure)
        return base
