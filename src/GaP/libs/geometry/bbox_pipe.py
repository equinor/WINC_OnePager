""" pipe builder
"""

import numpy as np

from ..models import (
    PipeCementModel,
)

from .bbox_utils import BBoxXY, BBoxZ

def bbox_for_pipe (cg: PipeCementModel,
                    LGR_sizes_xy: list[float], 
                    LGR_depths: np.ndarray, 
                    min_grd_size: float) -> tuple:
    """ build bounding box for the pipe

        Args:
            cg (PipeCementModel): contains geometrical information of the casing
            LGR_sizes_xy (list[float]): LGR xy grid intervals
            LGR_depths (np.ndarray): specific size of each LGR grid in the Z direction
            min_grd_size (float): minimum grid size

        Returns:

            tuple: the 3D bounding box and ID
    """

    # for convenience
    ID = cg.ID

    # z
    # pipe
    k_min_pipe, k_max_pipe = BBoxZ(cg.pipe.strt_depth, cg.pipe.end_depth, LGR_depths)
    
    # x-y
    x_min_pipe, x_max_pipe, y_min_pipe, y_max_pipe = BBoxXY(ID, LGR_sizes_xy, min_grd_size)

    return ID, x_min_pipe, x_max_pipe, y_min_pipe, y_max_pipe, k_min_pipe, k_max_pipe
