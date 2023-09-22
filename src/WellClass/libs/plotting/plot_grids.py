
from typing import Union

import numpy as np
import pandas as pd

from ..grid_utils.grid_coarse import GridCoarse
from ..grid_utils.grid_refine import GridRefine
from ..grid_utils.grid_lgr import GridLGR

from .plot_utils import plot_well_perm

def plot_grid(my_well, 
              grid: Union[GridCoarse, GridRefine, GridLGR], 
              *, 
              on_coarse):
    """ Plot well sketch and 2D slice of the permeability
    """

    # x-z corner coordinates at center y
    xcorn, zcorn = grid.extract_xz_corn_coords()

    # extract z PERM values at center y coarse grid
    Z = grid.extract_xz_slice()

    # plot x-z slice
    plot_well_perm(my_well, x=xcorn, y=zcorn, Z=Z, on_coarse=on_coarse)

def plot_coarse(my_well, grid_coarse: GridCoarse, *, on_coarse=True):
    """ Plot well sketch and 2D slice of the permeability, at coarse grid
    """

    # x-z corner coordinates at center y
    xcorn, zcorn = grid_coarse.extract_xz_corn_coords()

    # extract z PERM values at center y coarse grid
    Z = grid_coarse.extract_xz_slice()

    # plot x-z slice
    plot_well_perm(my_well, x=xcorn, y=zcorn, Z=Z, on_coarse=on_coarse)

def plot_refine(my_well, grid_refine: GridRefine):
    """ Visualization of well sketch and LGR grids
    """

    # x-z corner coordinates at center y
    xcorn, zcorn = grid_refine.extract_xz_corn_coords()

    # LGR
    # extract z PERM values at center y coarse grid
    Z = grid_refine.extract_xz_slice()

    # plot it
    plot_well_perm(my_well, x=xcorn, y=zcorn, Z=Z, on_coarse=False)  