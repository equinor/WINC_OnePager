""" pydantic model for pipe and cement
"""

# handle type hints problem for python version < 3.10
from typing import Union

from pydantic import BaseModel

from .model_depth import DepthModel

class ElemModel(BaseModel):
    """ define base model
    """
    ID: float
    pipe: DepthModel
    type: Union[str, None] = None
    
class PipeCementModel(ElemModel):

    cement: Union[DepthModel, None] = None
    oph: Union[DepthModel, None] = None
