
import pandas as pd

from src.GaP.libs.carfin import build_grdecl

from .grid_coarse import GridCoarse

from .df2gap import (
    to_gap_casing_list,
    to_gap_barrier_list
)

from .LGR_builder_base import LGRBuilderBase

class LGRBuilder(LGRBuilderBase):

    def __init__(self,
                 grid_init: GridCoarse, 
                 annulus_df: pd.DataFrame,
                 drilling_df: pd.DataFrame, 
                 casings_df: pd.DataFrame, 
                 borehole_df: pd.DataFrame,
                 barriers_mod_df: pd.DataFrame,
                 Ali_way: bool):
        """ LGR grid information in x, y, z directions. We are going to compute the grid sizes in lateral (x and y) and vertical directions

            Args:
                grid_init (GridCoarse): all information about coarse grid
                annulus_df (pd.DataFrame): information about annulus
                drilling_df (pd.DataFrame): information about drilling
                casings_df (pd.DataFrame): information about casings and cement-bond
                borehold_df (pd.DataFrame): information about borehole
                barriers_mod_df (pd.DataFrame): information about barrier
                Ali_way (bool): use Ali's algorithm to compute lateral grids and apply refdepth in z direction
        """

        super().__init__(
            grid_init,
            annulus_df,
            drilling_df,
            casings_df,
            borehole_df,
            barriers_mod_df,
            Ali_way
        )

    def build_grdecl(self, 
                     output_dir: str, output_name: str,
                     drilling_df: pd.DataFrame, 
                     casings_df: pd.DataFrame, 
                     barriers_mod_df: pd.DataFrame):
        """ build grdecl file and output it
        """

        # prepare info about Casing, Cement Bond and Open hole  for GaP
        casing_list = to_gap_casing_list(drilling_df, 
                                         casings_df)

        # prepare info about Barrier for GaP 
        barrier_list = to_gap_barrier_list(barriers_mod_df)

        # generate .grdecl file
        # TODO(hzh): add 1s to indices here
        build_grdecl(output_dir, output_name,
                     casing_list,
                     barrier_list,
                     self.LGR_sizes_x, 
                     self.LGR_depths, 
                     self.LGR_numb_z, 
                     self.min_grd_size,
                     self.NX, self.NY,
                     self.main_grd_i + 1, self.main_grd_j + 1,
                     self.main_grd_min_k + 1, self.main_grd_max_k + 1,
                     self.no_of_layers_in_OB)


