""" BoundingBox of element in depth direction
"""

from pydantic import BaseModel

class DepthModel(BaseModel):
    """ Contains the start/end depths and permeability value

        Args:
            strt_depth (float): starting depth of element (pipe/cement-bond/openhole)
            end_depth (float): ending depth of element
            perm (float|int): permeability, optional
    """
    strt_depth: float
    end_depth: float
    perm: float|int|None = None