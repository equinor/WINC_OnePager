from pydantic import BaseModel

from .well_models_utils import HoleCasingModel, PlugsModel, StratigraphyModel, SubsurfaceAssumptionsModel, WellHeaderModel, WellSurveyModel


class MetaDataModel(BaseModel):
    """
    meta data

    Args:
        namespace (str): name space
        name (str): can use it for project name
        author (str): who made this yaml file

    """

    namespace: str = "screen"
    name: str | None = None
    author: str | None = None


class WellSpec(BaseModel):
    """
    specs for standard well information
    Args:
        well_header (WellHeaderModel): well header information
        well_survey (WellSurveyModel): well survey information
        hole_casings (list[HoleCasingModel]): list of hole casing information
        plugs (list[PlugsModel]): list of plug information
        stratigraphy (list[StratigraphyModel]): list of stratigraphy information
        subsurface_assumptions (SubsurfaceAssumptionsModel): subsurface assumptions information
        co2_datum (CO2DatumModel): co2 datum
    """

    well_header: WellHeaderModel
    well_survey: WellSurveyModel | None = None
    hole_casings: list[HoleCasingModel] | None = None
    plugs: list[PlugsModel] | None = None
    stratigraphy: list[StratigraphyModel] | None = None


class WellModellingSpec(WellSpec):
    """
    extra specs for pressure information
    Args:
        reservoir_pressure (ReservoirPressureModel): general reservoir pressure information
        main_barrier (str): main barrier name to compute pressure
    """

    subsurface_assumptions: SubsurfaceAssumptionsModel | None = None


class WellModel(BaseModel):
    """
    model including all parameters
    Args:
        apiVersion (str): current version of this yaml format
        kind (str): for GaP
        metadata (MetaDataModel): miscelaneous data
        spec (WellPressureSpec): well specification
    """

    apiVersion: str = "well/v0.1"
    kind: str = "Well"
    metadata: MetaDataModel | None = None
    spec: WellModellingSpec
