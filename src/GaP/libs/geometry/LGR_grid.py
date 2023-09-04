""" generate fine grid from well casing geometries
"""

import numpy as np

from .LGR_grid_utils import (
    compute_min_grid_size,
    compute_max_num_of_fine_grid_xy,
    compute_LGR_xy,
    compute_LGR_z,
)

def generate_LGR_grids(casing_IDs: list[float],
                        main_grd_dx, main_DZ,
                        ref_depth: float,
                        main_grd_min_k: int, main_grd_max_k: int, 
                        no_of_layers_in_OB: int,
                        use_default_grid_size=True) -> tuple[list[float], np.ndarray, np.ndarray, float]:
    """ generate LGR grid information

    Args:

        casing_IDs (list[float]): list of IDs of conductor, surface and production
        main_grd_dx (float): dx for the cell (coarse grid) at the well location 
        main_DZ (np.ndarray): thickness of each layer in the Z direction (coarse grid)
        main_grd_min_k (int): minimum z index (coarse grid)
        main_grd_max_k (int): maximum z index (coarse grid)
        no_of_layers_in_OB (int): number of layers in Overburden
        use_default_grid_size (bool): flag for using default grid size, which is 5cm 

    Returns:

        LGR_sizes_xy (list[float]): LGR xy grid intervals
        LGR_depths (np.ndarray): specific size of each LGR grid in the Z direction
        LGR_numb_z (np.ndarray): number of chops for each main DZ
        min_grd_size (float): minimum dz value (LGR grid)
    """

    # for convenience
    cond_Casing_ID = casing_IDs[0]
    surf_Casing_ID = casing_IDs[1]
    prod_Casing_ID = casing_IDs[2]

    # 0. grid size for fine grid
    if use_default_grid_size:

        # TODO(hzh): fixed value???????
        min_grd_size = 0.05

    else:
        # compute min grid size
        min_grd_size = compute_min_grid_size(cond_Casing_ID, 
                                             surf_Casing_ID, 
                                             prod_Casing_ID)

    # 1. number of fine grids for the well core
    no_latral_fine_grd = compute_max_num_of_fine_grid_xy(cond_Casing_ID, 
                                                         surf_Casing_ID, 
                                                         prod_Casing_ID, 
                                                         min_grd_size)
    # 2. refine grid in x-y directions
    LGR_sizes_xy = compute_LGR_xy(no_latral_fine_grd, 
                                    main_grd_dx,
                                    min_grd_size)

    # 3. refine grid in z direction
    LGR_depths, LGR_numb_z = compute_LGR_z(main_DZ,
                                            ref_depth,
                                            main_grd_min_k, main_grd_max_k, 
                                            no_of_layers_in_OB)

    return LGR_sizes_xy, LGR_depths, LGR_numb_z, min_grd_size