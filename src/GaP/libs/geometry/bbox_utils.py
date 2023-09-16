""" compute bounding boxes for lateral x-y and depth z directions
"""
from typing import Tuple, List

import math
import numpy as np

def BBoxXY(ID: float,
            LGR_sizes_xy: List[float], 
            min_grd_size: float
           ) -> Tuple[int, int, int, int]:
    """ Compute bounding box in x-y directions

        Args:
            ID (float): Internal diameter value
            LGR_sizes_xy (list[float]): LGR xy grid intervals
            min_grd_size (float): minimum grid size

        Returns:

            tuple: the 2D bounding box, indices, in lateral x-y plane of LGR grid
    """

    # x
    no_grd_elem_x = math.floor(ID/min_grd_size)

    x_min_elem = math.ceil((len (LGR_sizes_xy)-no_grd_elem_x)/2)
    if x_min_elem < 0:
        raise ValueError('ERROR:The ID of barrier is larger than refined area')
    x_max_elem = x_min_elem + no_grd_elem_x 
    
    # y
    no_grd_elem_y = math.floor(ID/min_grd_size)

    y_min_elem = math.ceil((len (LGR_sizes_xy)-no_grd_elem_y)/2)
    if y_min_elem < 0:
        raise ValueError('ERROR:The ID of tube is larger than refined area')
    y_max_elem = y_min_elem + no_grd_elem_y 

    return x_min_elem, x_max_elem, y_min_elem, y_max_elem

def BBoxZ(strt_depth: float,
          end_depth: float,
          LGR_depths: np.ndarray
          ) -> Tuple[int, int]:
    """ Compute bounding box in z direction

        Args:
            strt_depth (float): start depth of element, such as pipe, cement-bond, or barrier
            end_depth (float): end depth of element
            LGR_depths (np.ndarray): specific size of each LGR grid in the Z direction

        Returns:
            tuple: min/max z indices in LGR grid
    """

    # pipe
    if (LGR_depths[0] - strt_depth) > 0:
        k_min_elem = 1 
    else:
        k_min_elem = np.argmin(abs(LGR_depths - strt_depth))+ 1 

    if (LGR_depths[-1] - end_depth) < 0:
        k_max_elem =  len(LGR_depths)
    else:
        k_max_elem = np.argmin(abs(LGR_depths - end_depth))+ 1 

    return k_min_elem, k_max_elem
