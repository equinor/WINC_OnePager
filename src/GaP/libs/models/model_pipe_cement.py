""" pydantic model for pipe and cement
"""

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
