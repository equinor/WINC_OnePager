
# handle type hints problem for python version < 3.10
from typing import Union, List

from pydantic import BaseModel

from .well_model_utils import (
    WellHeaderModel,
    DrillingModel,
    CasingCementModel,
    BarrierModel,
    BarrierPermeabilityModel,
    GeologyModel,
    AssumptionsModel,
    ReservoirPressureModel,
    # CO2DatumModel,
    # MainBarrierModel,
)

class MetaDataModel(BaseModel):
    """ meta data
        Args:
            namespace (str): name space
            name (str): can use it for project name
            author (str): who made this yaml file
    """
    namespace: str = 'screen'
    name: Union[str, None] = None
    author: Union[str, None] = None

class WellSpec(BaseModel):
    """ specs for standard well information
        Args:
            well_header (WellHeaderModel): well header information
            drilling (list[DrillingModel]): list of drilling information
            casing_cement (list[CasingCementModel]): list of casing information
            barrier (list[BarrierModel]):  list of barrier information
            barrier_permeability (BarrierPermeabilityModel): list of barrier permeability
            geology (list[GeologyModel]): list of geology, such as formations, information
            assumptions (AssumptionsModel): misceleaneous information
            co2_datum (CO2DatumModel): co2 datum 
    """
    well_header: WellHeaderModel
    drilling: List[DrillingModel]
    casing_cement: List[CasingCementModel]
    barriers: Union[List[BarrierModel], None] = None
    barrier_permeability: Union[BarrierPermeabilityModel, None] = None
    geology: Union[List[GeologyModel], None] = None
    assumptions: Union[AssumptionsModel, None] = None
    co2_datum: Union[float, int]                                       # CO2DatumModel

class WellPressureSpec(WellSpec):
    """ extra specs for pressure information
        Args:
            reservoir_pressure (ReservoirPressureModel): general reservoir pressure information
            main_barrier (str): main barrier name to compute pressure
    """
    reservoir_pressure: Union[ReservoirPressureModel, None] = None
    main_barrier: Union[str, None] = None

class WellModel(BaseModel):
    """ model including all parameters
        Args:
            apiVersion (str): current version of this yaml format
            kind (str): for GaP
            metadata (MetaDataModel): miscelaneous data
            spec (WellPressureSpec): well specification
    """
    apiVersion: str = 'well/v0.1'
    kind: str = 'Well'
    metadata: Union[MetaDataModel, None] = None
    spec: WellPressureSpec
