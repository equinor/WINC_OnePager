""" barrier builder 
"""

import numpy as np

from ..models import (
    ElemModel,
)

from .bbox_utils import BBoxXY, BBoxZ

def bbox_for_barrier (cg: ElemModel, 
                      LGR_sizes_xy: list[float], 
                      LGR_depths: np.ndarray, 
                      min_grd_size: float) -> tuple:
    """ build bouding box for barrier

        Args:
            cg (ElemModel): contains geometrical information of the casing
            LGR_sizes_xy (list[float]): LGR xy grid intervals
            LGR_depths (np.ndarray): specific size of each LGR grid in the Z direction
            min_grd_size (float): minimum grid size

        Returns:

            tuple: the 3D bounding box and ID

    """

    # for convenience
    ID = cg.ID

    # z
    k_min_bar, k_max_bar = BBoxZ(cg.pipe.strt_depth, cg.pipe.end_depth, LGR_depths)
    
    # x-y
    x_min_bar, x_max_bar, y_min_bar, y_max_bar = BBoxXY(ID, LGR_sizes_xy, min_grd_size)

    return ID, x_min_bar, x_max_bar, y_min_bar, y_max_bar, k_min_bar, k_max_bar

