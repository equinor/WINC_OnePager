
from typing import Tuple

import numpy as np

# 0. Compute number of cells from the diameter and min grid size
def compute_ngrd(diam: float, min_grd_size: float) -> int:
    """ Computes floor and ceiling of discretizing diameter

    Args:
        diam (float): diameter
        min_grd_size (float): minimum grid size

    Returns:
        n_grd (int): the discretization with smallest discrepancy
    """
    
    area = np.pi * (diam/2)**2
    cart_radius = np.sqrt(area) / 2

    ceil_nx =  2*int(np.ceil(cart_radius/ min_grd_size))
    floor_nx = 2*int(np.floor(cart_radius / min_grd_size))
                
    ceil_diff = np.abs(area - (ceil_nx*min_grd_size)**2)
    floor_diff = np.abs(area -(floor_nx*min_grd_size)**2)

    if floor_diff < ceil_diff:
        n_grd = floor_nx
    else:
        n_grd = ceil_nx

    return n_grd

# 1. Compute LGR grids in x-y directions
def generate_LGR_xy(no_latral_fine_grd: int, 
                    min_grd_size: float, 
                    main_grd_dx: float, 
                    main_grd_dy: float,
                    *,
                    Ali_way: bool=False):
    """ generate LGR grid in x-y directions

        Args:
            no_latral_fine_grd (int): number of lateral grids (fine grid)
            min_grd_size (float): minimum grid size, default 5cm
            main_grd_dx (float): x grid size (coarse grid)
            main_grd_dy (float): y grid size (coarse grid)
            Ali_way (bool): use Ali's algorithm

        Returns:
            tuples: LGR_sizes_x, LGR_sizes_y, DX_log_inc
    """
    
    if Ali_way:

        # add two more at both end
        no_latral_fine_grd += 4

        # add transistion zone between coarse and fine grids
        LGR_size_no_outer = [min_grd_size*100] + \
                            [min_grd_size*10] + \
                            [min_grd_size] * no_latral_fine_grd + \
                            [min_grd_size*10] + \
                            [min_grd_size*100]

        # 6. LGR_sizes_xy, adding dx on both sides to make it one cell size
        LGR_sizes_xy = [(main_grd_dx - sum(LGR_size_no_outer))/2] + LGR_size_no_outer + [(main_grd_dx - sum(LGR_size_no_outer))/2]

        return LGR_sizes_xy, LGR_sizes_xy, None
    
    else:
        # 1.1 generate a LGR array by repeating min_grd_size with the number of fine cells
        LGR_size_fine_grd = [min_grd_size]*no_latral_fine_grd

        # 1.2. compute transition zone between finer and coarse meshes
        #Number of logarithmic increments between smallest mesh and coarse mesh
        n_log_inc = 3
        
        #Distance between cells representing the wells and coarse cells
        DX_transition = (main_grd_dx - (no_latral_fine_grd*min_grd_size)) / 2
        DY_transition = (main_grd_dy - (no_latral_fine_grd*min_grd_size)) / 2
        
        #Compute corner edges of logarithmic increments
        Xcorn_log_inc = np.logspace(np.log10(min_grd_size), np.log10(DX_transition + min_grd_size), num=n_log_inc+1) 
        Ycorn_log_inc = np.logspace(np.log10(min_grd_size), np.log10(DY_transition + min_grd_size), num=n_log_inc+1) 
        
        #Compute DX DY of log increments
        DX_log_inc = np.diff(Xcorn_log_inc)
        DY_log_inc = np.diff(Ycorn_log_inc)
        
        # 1.3. combines transition zones and inner zone around well
        LGR_sizes_x = np.concatenate((DX_log_inc[::-1], LGR_size_fine_grd, DX_log_inc))
        LGR_sizes_y = np.concatenate((DY_log_inc[::-1], LGR_size_fine_grd, DY_log_inc))

        return LGR_sizes_x, LGR_sizes_y, DX_log_inc

# 2. compute LGR grids in z direction
def generate_LGR_z(DZ_rsrv: np.ndarray, 
                   DZ_ovb_coarse: np.ndarray,
                   ref_depth: float=0.0, 
                   z_finer_scale=10) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """ generate LGR grid in z direction

        Args:
            DZ_rsrv (np.ndarray): dz values for reservoir
            DZ_ovb_coarse (np.ndarray): dz values for overburden
            ref_depth (float): depth where LGR Z starts, default: 0
            z_finer_scale (int): LGR for ovb on finer grid, default: 10, by making it 10 times smaller

        Returns:
            tuples: LGR_sizes_z, LGR_numb_z, LGR_depths, DZ_ovb
    """
    
    # 2.2 for ovb

    # refine grid
    dz_finer = DZ_ovb_coarse/z_finer_scale
    DZ_ovb = np.repeat(dz_finer, z_finer_scale)
    
    # 2.3 combine ovb and rsrv
    LGR_sizes_z = np.concatenate((DZ_ovb, DZ_rsrv))

    # ratio
    LGR_numb_z = len(DZ_ovb_coarse)*[z_finer_scale] + len(DZ_rsrv)*[1]
    LGR_numb_z = np.array(LGR_numb_z)

    # add ref_depth
    LGR_sizes_z[0] += ref_depth

    # accumulation depths
    LGR_depths = LGR_sizes_z.cumsum()

    return LGR_sizes_z, LGR_numb_z, LGR_depths, DZ_ovb
