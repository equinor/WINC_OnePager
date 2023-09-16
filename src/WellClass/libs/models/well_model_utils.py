
# handle type hints problem for python version < 3.10
from typing import Union, List

from pydantic import BaseModel, validator

from ..utils.fraction_float import fraction_float

class WellHeaderModel(BaseModel):
    """ General information about the well.
    
        Args:
            well_name (str): well name
            well_rkb (float): RKB elevation
            sf_depth_msl (float): depth of the sea floor
            well_td_rkb (float): total depth of RKB
            sf_temp (float): sea floor temperature
            geo_tgrad  (float): geothermal gradient
    """
    well_name: str
    well_rkb: Union[int, float]
    sf_depth_msl: Union[int, float]
    well_td_rkb: Union[int, float]
    sf_temp: Union[int, float]
    geo_tgrad: Union[int, float]

class DrillingRawModel(BaseModel):
    """ Information about the drilling intervals of the well

        Args:
            top_rkb (float): the top depth in RKB
            bottom_rkb (float): the bottom depth in RKB 
            diameter_in (float, str): the diameter of the borehole in inches
    """
    top_rkb: Union[int, float]
    bottom_rkb: Union[int, float]
    diameter_in: Union[float, int, str]

    @validator('diameter_in')
    def diameter_in_converter(cls, v):
        if isinstance(v, (float, int)):
            return v
        elif isinstance(v, str):
            return fraction_float(v)
        else:
            raise ValueError('diameter_in must be a float or string')
        
class DrillingModel(DrillingRawModel):
    """ Information about the drilling intervals of the well

        Args:
            oh_perm (float): faked permeability for open-hole 
    """
    oh_perm: Union[int, float] = 10000


class CasingCementModel(DrillingRawModel):
    """ Information about the casing and cementing intervals of the well

        Args:
          toc_rkb (float): top of cement-bond in RKB  
          boc_rkb (float): bottom of cement-bond in RKB 
          shoe (bool): whether or not it has a shoe
          cb_perm (float): permeability for cement-bond
    """
    toc_rkb: Union[int, float]
    boc_rkb: Union[int, float]
    shoe: bool
    cb_perm: Union[int, float] = 5

class BarrierModel(BaseModel):
    """ Information about the barrier in the well 

        Args:
            barrier_name (str): the barrier name 
            barrier_type (str): the barrier type 
            top_rkb (float): the top depth in RKB
            bottom_rkb (float): the bottom depth in RKB
            barrier_perm (float): permeability for the barrier
    """
    barrier_name: str
    barrier_type: str
    top_rkb: Union[int, float]
    bottom_rkb: Union[int, float]
    barrier_perm: Union[int, float] = 0.5

class GeologyModel(BaseModel):
    """ The geological units encountered in the well

        Args:
            top_rkb (float): top depth in RKB
            geol_unit (str): name 
            reservoir_flag (bool): whether or not it is considered a reservoir  
    """
    top_rkb: Union[int, float]
    geol_unit: str
    reservoir_flag: bool

class AssumptionsModel(BaseModel):
    """ model for assumptions
    """
    pass

class ReservoirPressureModel(BaseModel):
    """ Model for reservoir pressure

        Args:
            depth_msl (float): mean sea level depth
            RP1 (str): reservoir pressure
            RP2 (str): reservoir pressure
            RP3 (str): reservoir pressure
    """
    depth_msl: Union[int, float]
    RP1: Union[str, float, int, None] = None
    RP2: Union[str, float, int, None] = None
    RP3: Union[str, float, int, None] = None

    # @validator('RP1', 'RP2', 'RP3')
    # def diameter_in_converter(cls, v):
    #     if isinstance(v, float, int):
    #         return v
    #     elif isinstance(v, str):
    #         return float(v)
    #     else:
    #         raise ValueError('diameter_in must be a float or string')
        
class CO2DatumModel(BaseModel):
    """ Model for CO2 datum. Note it is not used.

        Args:
            co2_msl (float): CO2 datum depth
    """
    co2_msl: Union[int, float]

class MainBarrierModel(BaseModel):
    """ model for main barrier

        Args:
            barrier_name (str): barrier name for proxy compute
    """
    barrier_name: str

class BarrierPermeabilityModel(BaseModel):
    """ model for Barrier permeability

        Args:
            quality (str): list of quality level
            kv (float): list of permeability values
    """
    quality: Union[List[str], None] = None
    kv: List[Union[float, int]]

