""" Module to output the .grdecl file for LGR grid
"""
from typing import List

import os
import sys

import numpy as np

from ..models import (
    PipeCementModel,
    ElemModel
)

from .CARFIN_core import (
    pre_CARFIN,
    CARFIN_keywords,
    coreCARFIN,
    endCARFIN
)

def build_grdecl(output_folder: str, 
                 LGR_NAME: str,
                 casing_list: List[PipeCementModel],
                 barrier_list: List[ElemModel], 
                 LGR_sizes_xy: List[float], 
                 LGR_depths: np.ndarray,
                 LGR_numb_z: np.ndarray,
                 min_grd_size: float,
                 NX: int, NY: int,
                 main_grd_i: int, main_grd_j: int,
                 main_grd_min_k: int, main_grd_max_k: int,
                 no_of_layers_in_OB: int) -> None:
    """ Generate grdecl grid

        To be able to model a legacy well in reservoir scale, we need to make sure all of the elements including\
        multiple casings with different OD, cement bonds, barriers, the open area between barriers inside the casings, etc.\
        are considered. 

        Then, the existing simple well models available in the commercial simulators (Eclipse, PFT, IX) are not  able to\ 
        include those details. They just introduce a node where the flow will be discharged from/to the grid to the node. 

        The way around it is to define the well as a part of the reservoir by manipulating the local grids (LGR) and properties of the grid, 

        Args:
            output_folder (str): output directory
            LGR_NAME (str): output filename, without the suffix '.grdecl'
            casing_list (list[PipeCementModel]): contains list of casing geometry
            barrier_list (list[ElemModel]): contains list of barriers
            LGR_sizes_xy (list[float]): LGR xy grid intervals
            LGR_depths (np.ndarray): specific size of each LGR grid in the Z direction
            LGR_numb_z (np.ndarray): number of chops for each main DZ
            min_grd_size (float): minimum grid size
            NX (int): number of grids in x direction (coarse grid)
            NY (int): number of grids in y direction (coarse grid)
            main_grd_i (int): x index of well location in coarse grid
            main_grd_j (int): y index of well location in coarse grid
            main_grd_min_k (int): minimum z index in coarse grid
            main_grd_max_k (int): maximum z index in coarse grid
            no_of_layers_in_OB (int): number of layers in overburden (coarse grid)

        Returns:
            None             
    """

    # 0. prepare file for output

    # check output directory
    if not os.path.exists(output_folder):
        os.makedirs(output_folder, exist_ok=True)

    # generate output file name
    out_fname = os.path.join(output_folder, LGR_NAME+'.grdecl')

    # open it
    if os.path.exists(out_fname) == True:
        O = open(out_fname,"r+")
    else: 
        O = open(out_fname,"x")

    O.truncate(0)

    # 1. start the process
    pre_CARFIN(LGR_NAME,
                NX, NY,
                main_grd_i, main_grd_j, 
                no_of_layers_in_OB, 
                O)
    
    CARFIN_keywords(LGR_NAME,
                    main_grd_i, main_grd_j, 
                    main_grd_min_k, main_grd_max_k, 
                    LGR_sizes_xy, 
                    LGR_numb_z, 
                    min_grd_size,
                    O)

    coreCARFIN(LGR_NAME,
                casing_list,
                barrier_list,
                LGR_sizes_xy, 
                LGR_depths, 
                min_grd_size,
                O)

    # find minimum IDs
    reopen_ID = sys.float_info.max
    for casing in casing_list:
        if casing.ID < reopen_ID:
            reopen_ID = casing.ID

    endCARFIN(LGR_NAME,
                reopen_ID,
                LGR_sizes_xy, 
                main_grd_min_k, 
                min_grd_size,
                no_of_layers_in_OB,
                O)

    # 2. done
    O.close()

    # for qc
    print ('Output LGR CARFIN to: ', os.path.abspath(out_fname))