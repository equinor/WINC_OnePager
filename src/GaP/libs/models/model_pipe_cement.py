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

class PipeModel(ElemModel):
    """ define model for pipe
    """
    oph: Union[DepthModel, None] = None

class PipeCementModel(PipeModel):
    """ define model for pipe and cement
    """
    cement: Union[DepthModel, None] = None

