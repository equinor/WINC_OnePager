
from typing import Tuple

import numpy as np
import pandas as pd

from .grid_coarse import GridCoarse
from .grid_refine_base import GridRefineBase

from .extract_grid_utils import (
    extract_xz_corn_coords,
    extract_xz_prop_slice,
)

class GridRefine(GridRefineBase):

    def __init__(self, 
                 grid_coarse: GridCoarse,
                 LGR_sizes_x, LGR_sizes_y, LGR_sizes_z,
                 min_grd_size: float
                 ):
        """ dataframe for LGR mesh for the center coarse grid
        """

        super().__init__(grid_coarse, 
                         LGR_sizes_x, LGR_sizes_y, LGR_sizes_z,
                         min_grd_size)
        
    def build_LGR(self, 
                  drilling_df: pd.DataFrame, 
                  casings_df: pd.DataFrame, 
                  barriers_mod_df: pd.DataFrame) -> None:
        """ assign material types to corresponding permeabilities

            Args:

                drilling_df (pd.DataFrame): information about drilling
                casings_df (pd.DataFrame): information about casings and cement-bond
                borehold_df (pd.DataFrame): information about borehole
                barriers_mod_df (pd.DataFrame): information about barrier                  
        """

        # 1. compute lateral number of refined grid
        self._compute_num_lateral_fine_grd(drilling_df, casings_df, barriers_mod_df)

        # 2. set bounding box
        self._compute_bbox(drilling_df, casings_df, barriers_mod_df)

        # 3. set material type
        self._set_material_type(drilling_df, casings_df, barriers_mod_df)

        # 4. set permeability
        self._set_permeability(drilling_df, casings_df, barriers_mod_df)

 
    def extract_xz_corn_coords(self) -> Tuple[np.ndarray, np.ndarray]:
        """ generate xcorn and zcorn coordinates
        """

        # for convenience
        mesh_df = self.mesh_df

        # for shifting
        sDX = self.main_grd_dx/2
        sDY = self.main_grd_dy/2

        # generate grid coordinates for plotting
        xcorn, zcorn = extract_xz_corn_coords(mesh_df, sDX, sDY)

        return xcorn, zcorn
    
    def extract_xz_slice(self) -> np.ndarray:
        """ generate x-z PERM slice
        """
        # for convenience
        mesh_df = self.mesh_df

        # extract permeability
        Z = extract_xz_prop_slice(mesh_df)

        return Z 
