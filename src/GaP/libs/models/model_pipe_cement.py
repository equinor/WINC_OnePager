""" pydantic model for pipe and cement
"""

# handle type hints problem for python version < 3.10
from __future__ import annotations

from pydantic import BaseModel

from .model_depth import DepthModel

class ElemModel(BaseModel):
    """ define base model
    """
    ID: float
    pipe: DepthModel
    type: str|None = None
    
class PipeCementModel(ElemModel):

    cement: DepthModel|None = None
    oph: DepthModel|None = None
