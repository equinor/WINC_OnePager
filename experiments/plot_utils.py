
import numpy as np
import pandas as pd

from src.WellClass.libs.plotting.plot_utils import plot_well_perm

from .grid_coarse import GridCoarse
from .grid_refine import GridRefine

def plot_coarse(my_well, grid_coarse: GridCoarse):
    """ Plot well sketch and 2D slice of the permeability, at coarse grid
    """

    # x-y corner coordinates
    xcorn, ycorn = grid_coarse.extract_xy_corn_coords()

    # extract z PERM values at center y coarse grid
    Z = grid_coarse.extract_xz_slice()

    # plot x-z slice
    plot_well_perm(my_well, x=xcorn, y=ycorn, Z=Z, on_coarse=True)

def plot_refine(my_well, grid_refine: GridRefine):
    """ Visualization of well sketch and LGR grids
    """

    # x-y corner coordinates
    xcorn, ycorn = grid_refine.xtract_xy_corn_coords()

    # LGR
    # extract z PERM values at center y coarse grid
    Z = grid_refine.extract_xz_slice()

    # plot it
    plot_well_perm(my_well, x=xcorn, y=ycorn, Z=Z, on_coarse=False)  