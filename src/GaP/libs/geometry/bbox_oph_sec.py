""" pipe builder with open hole sections
"""
from typing import List

import numpy as np

from ..models import (
    ElemModel,
)

from .bbox_utils import BBoxXY, BBoxZ

def bbox_for_oph_sec (cg: ElemModel,
                        LGR_sizes_xy: List[float], 
                        LGR_depths: np.ndarray, 
                        min_grd_size: float) -> tuple:
    """ build bounding box for pipe with open hole sections

        Args:
            cg (ElemModel): contains geometrical information of the open-hole section
            LGR_sizes_xy (list[float]): LGR xy grid intervals
            LGR_depths (np.ndarray): specific size of each LGR grid in the Z direction
            min_grd_size (float): minimum grid size

        Returns:

            tuple: the 3D bounding box and min/max k hole and ID
    """

    # for convenience
    ID = cg.ID

    # z, open hole section
    k_min_hole, k_max_hole = BBoxZ(cg.pipe.strt_depth, cg.pipe.end_depth, LGR_depths)
    
    # x-y
    x_min_pipe, x_max_pipe, y_min_pipe, y_max_pipe = BBoxXY(ID, LGR_sizes_xy, min_grd_size)
    
    return ID, x_min_pipe, x_max_pipe, y_min_pipe, y_max_pipe, k_min_hole, k_max_hole

