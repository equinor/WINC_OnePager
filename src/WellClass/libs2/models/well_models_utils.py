# handle type hints problem for python version < 3.10
from typing import Literal

from pydantic import BaseModel, field_validator

from ..utils.fraction_float import fraction_float


class WellHeaderModel(BaseModel):
    """
    General information about the well.

    Args:
        well_name (str): well name
        well_rkb (float): RKB elevation
        sf_depth_msl (float): depth of the sea floor
        well_td_rkb (float): total depth of RKB
        sf_temp (float): sea floor temperature
        geo_tgrad (float): geothermal gradient

    """

    unique_wellbore_identifier: str
    depth_reference_rkb: int | float
    depth_reference_rkb_unit: Literal["ft", "m"]
    ground_elevation: int | float
    ground_elevation_unit: Literal["ft", "m"]
    total_depth_rkb: int | float
    total_depth_rkb_unit: Literal["ft", "m"]


class WellSurveyModel(BaseModel):
    """
    Information about the well survey

    Args:
        md_rkb (float): measured depth in RKB
        inclination_deg (float): inclination in degree
        azimuth_deg (float): azimuth in degree

    """

    md_rkb: list[int | float]
    inclination_deg: list[int | float]
    azimuth_deg: list[int | float]


class HoleCasingModelRaw(BaseModel):
    """
    Information about the drilling intervals of the well

    Args:
        top_rkb (float): the top depth in RKB
        bottom_rkb (float): the bottom depth in RKB
        diameter_in (float, str): the diameter of the borehole in inches

    """

    name: str
    type: Literal["hole", "casing", "casing cement"]
    top_rkb: int | float
    bottom_rkb: int | float
    diameter_in: float | int | str
    shoe: bool | None = False

    @field_validator("diameter_in")
    def diameter_in_converter(cls, v: float | str) -> float:
        if isinstance(v, float):
            return v
        if isinstance(v, str):
            return fraction_float(v)
        raise ValueError("diameter_in must be a float or string")


class HoleCasingModel(HoleCasingModelRaw):
    """
    Information about the drilling intervals of the well

    Args:
        hc_perm (float): faked permeability for hole casing element

    """

    hc_perm: int | float | None = None


class PlugsRaw(BaseModel):
    """
    Information about the barrier in the well

    Args:
        barrier_name (str): the barrier name
        barrier_type (str): the barrier type
        top_rkb (float): the top depth in RKB
        bottom_rkb (float): the bottom depth in RKB
        barrier_perm (float): permeability for the barrier

    """

    name: str
    type: Literal["cement", "mechanical plug"]
    top_rkb: int | float
    bottom_rkb: int | float


class PlugsModel(PlugsRaw):
    """
    Information about the barrier in the well

    Args:
        barrier_perm (float): permeability for the barrier

    """

    cement_perm: int | float | None = None


class StratigraphyRaw(BaseModel):
    """
    The geological units encountered in the well

    Args:
        top_rkb (float): top depth in RKB
        geol_unit (str): name
        reservoir_flag (bool): whether or not it is considered a reservoir

    """

    name: str
    top_rkb: int | float
    bottom_rkb: int | float


class StratigraphyModel(StratigraphyRaw):
    """
    The geological units encountered in the well

    Args:
        reservoir_flag (bool): whether or not it is considered a reservoir

    """

    unit_type: str | None = None
    unit_perm: int | float | None = None


class ShminDataPoint(BaseModel):
    depth: float  # depth (e.g., meters)
    shmin: float  # minimum horizontal stress (e.g., MPa)


class SubsurfaceAssumptionsScenario(BaseModel):
    """model for subsurface assumptions"""

    temperature_gradient: int | float | None = None
    ground_temperature: int | float | None = None
    sg_brine: float | None = None
    sg_fluid: float | None = None
    fluid_type: str | None = None
    shmin_gradient: float | None = None
    shmin_data: list[ShminDataPoint] | None = None
    salinity: int | float | None = None
    z_fluid_contact: int | float | None = None
    p_fluid_contact: int | float | None = None
    z_resrv: int | float | None = None
    p_resrv: int | float | None = None
    z_msad: int | float | None = None
    p_delta: int | float | None = None


class SubsurfaceAssumptionsModel(BaseModel):
    """model for subsurface assumptions"""

    scenarios: list[SubsurfaceAssumptionsScenario] | None = None

    @field_validator("scenarios", mode="before")
    def ensure_list(cls, v: list[SubsurfaceAssumptionsScenario] | SubsurfaceAssumptionsScenario | None) -> list[SubsurfaceAssumptionsScenario] | None:
        if v is None:
            return v
        if isinstance(v, list):
            return v
        return [v]
