""" combine pipe withe cement bond
"""
from typing import List

import numpy as np
from typing import TextIO

from ..models import (
    PipeCementModel,
    ElemModel,
)

from ..geometry.bbox_pipe_with_oph_sec import bbox_for_pipe_with_oph_sec
from ..geometry.bbox_oph_sec import bbox_for_oph_sec
from ..geometry.bbox_cement_bond import bbox_for_cement_bond
from ..geometry.bbox_barrier import bbox_for_barrier

from .CARFIN_pipe_with_oph import CARFIN_pipe_with_oph
from .CARFIN_oph import CARFIN_oph
from .CARFIN_cement_bond import CARFIN_cement_bond
from .CARFIN_barrier import CARFIN_barrier

def CARFIN_pipe_and_cement_bond_builder(casing_geom: PipeCementModel, 
                                        LGR_sizes_xy: List[float], 
                                        LGR_depths: np.ndarray, 
                                        min_grd_size: float,
                                        LGR_NAME: str, 
                                        O: TextIO):  # noqa: E741
    """ pipe with cement bond

        Args:
            casing_geom (PipeCementModel): contains casing geometry information
            LGR_sizes_xy (list[float]): LGR xy grid intervals
            LGR_depths (np.ndarray): specific size of each LGR grid in the Z direction
            min_grd_size (float): minimum grid size
            LGR_NAME (str): output filename, without the suffix '.grdecl'
            O (TextIO): opened file handle            
    """
    # extract permeability
    pipe_perm = casing_geom.pipe.perm
    cement_perm = casing_geom.cement.perm

    # 1. for pipe
    # 1.1 bounding box for pipe
    ID, x_min_pipe, x_max_pipe, y_min_pipe, y_max_pipe, k_min_pipe, k_max_pipe, k_min_hole, k_max_hole = \
    bbox_for_pipe_with_oph_sec (casing_geom, 
                                LGR_sizes_xy, LGR_depths, 
                                min_grd_size)
    
    # 1.1 CARFIN output
    CARFIN_pipe_with_oph(ID,
                            x_min_pipe, x_max_pipe,
                            y_min_pipe, y_max_pipe,
                            k_min_pipe, k_max_pipe,
                            k_min_hole, k_max_hole,
                            pipe_perm,
                            LGR_NAME,
                            O)
    
    # 2. for cement bond
    # 2.1 bounding box for cement bond
    ID, x_min_pipe, x_max_pipe, y_min_pipe, y_max_pipe, k_min_CB, k_max_CB = \
    bbox_for_cement_bond (casing_geom, 
                            LGR_sizes_xy, LGR_depths, 
                            min_grd_size)
    # 2.2 CARFN output
    CARFIN_cement_bond(ID, 
                        x_min_pipe, x_max_pipe,
                        y_min_pipe, y_max_pipe,
                        k_min_CB, k_max_CB,
                        perm=cement_perm,
                        LGR_NAME=LGR_NAME, 
                        O=O)

def CARFIN_oph_builder(casing_geom: ElemModel, 
                        LGR_sizes_xy: List[float], 
                        LGR_depths: np.ndarray, 
                        min_grd_size: float,
                        LGR_NAME: str, 
                        O: TextIO):  # noqa: E741
    """ builder for open-hole section

        Args:
            casing_geom (ElemModel): contains open-hole section geometry information
            LGR_sizes_xy (list[float]): LGR xy grid intervals
            LGR_depths (np.ndarray): specific size of each LGR grid in the Z direction
            min_grd_size (float): minimum grid size
            LGR_NAME (str): output filename, without the suffix '.grdecl'
            O (TextIO): opened file handle            
    """
    # extract permeability
    pipe_perm = casing_geom.pipe.perm

    # 1. for pipe
    # 1.1 bounding box for pipe
    ID, x_min_pipe, x_max_pipe, y_min_pipe, y_max_pipe, k_min_hole, k_max_hole = \
    bbox_for_oph_sec(casing_geom, 
                     LGR_sizes_xy, LGR_depths, 
                     min_grd_size)
    
    # 1.1 CARFIN output
    CARFIN_oph(ID,
                x_min_pipe, x_max_pipe,
                y_min_pipe, y_max_pipe,
                k_min_hole, k_max_hole,
                pipe_perm,
                LGR_NAME,
                O)
        
def CARFIN_barrier_builder(barrier_geom: ElemModel, 
                            LGR_sizes_xy: List[float],
                            LGR_depths: np.ndarray,
                            min_grd_size: float,
                            LGR_NAME: str, 
                            O: TextIO):  # noqa: E741
    """ pipe with cement bond

        Args:
            barrier_geom (ElemModel): contains barrier geometry information
            LGR_sizes_xy (list[float]): LGR xy grid intervals
            LGR_depths (np.ndarray): specific size of each LGR grid in the Z direction
            min_grd_size (float): minimum grid size
            LGR_NAME (str): output filename, without the suffix '.grdecl'
            O (TextIO): opened file handle         
    """
    
    # extract permeability for barrier
    barrier_perm = barrier_geom.pipe.perm
    
    # 1. for pipe
    # 1.1 bounding box for pipe
    ID, x_min_bar, x_max_bar, y_min_bar, y_max_bar, k_min_bar, k_max_bar = \
    bbox_for_barrier (barrier_geom, 
                        LGR_sizes_xy, LGR_depths, 
                        min_grd_size)
    
    # 1.1 CARFIN output
    CARFIN_barrier(ID, 
                     x_min_bar, x_max_bar,
                     y_min_bar, y_max_bar,
                     k_min_bar, k_max_bar,
                     barrier_perm,
                     LGR_NAME, 
                     O)
