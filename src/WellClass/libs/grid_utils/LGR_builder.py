
import os
import pandas as pd

from src.GaP.libs.carfin.CARFIN_core import (
    pre_CARFIN,
    CARFIN_keywords,
    endCARFIN
)

from .grid_coarse import GridCoarse

from .df2gap_core import (
    df_to_gap_casing,
    df_to_gap_barrier
)

from .LGR_builder_base import LGRBuilderBase

class LGRBuilder(LGRBuilderBase):

    def __init__(self,
                 grid_init: GridCoarse, 
                 annulus_df: pd.DataFrame,
                 drilling_df: pd.DataFrame, 
                 Ali_way: bool):
        """ LGR grid information in x, y, z directions. We are going to compute the grid sizes in lateral (x and y) and vertical directions

            Args:
                grid_init (GridCoarse): all information about coarse grid
                annulus_df (pd.DataFrame): information about annulus
                drilling_df (pd.DataFrame): information about drilling
                Ali_way (bool): use Ali's algorithm to compute lateral grids and apply refdepth in z direction
        """

        super().__init__(
            grid_init,
            annulus_df,
            drilling_df,
            Ali_way
        )

    def build_grdecl(self, 
                     output_folder: str, LGR_NAME: str,
                     drilling_df: pd.DataFrame, 
                     casings_df: pd.DataFrame, 
                     barriers_mod_df: pd.DataFrame):
        """ build grdecl file and output it
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
                    self.NX, self.NY,
                    self.main_grd_i+1, self.main_grd_j+1, 
                    self.no_of_layers_in_OB, 
                    O)
        
        # 2. keywods
        CARFIN_keywords(LGR_NAME,
                        self.main_grd_i+1, self.main_grd_j+1, 
                        self.main_grd_min_k+1, self.main_grd_max_k+1, 
                        self.LGR_sizes_x, 
                        self.LGR_numb_z, 
                        self.min_grd_size,
                        O)

        # 3. the pipes/openholes/barriers
        df_to_gap_casing(drilling_df, 
                         casings_df,
                         LGR_NAME,
                         O)
        
        df_to_gap_barrier(barriers_mod_df,
                            LGR_NAME,
                            O)
        
        # find minimum IDs for casings
        reopen_ID = casings_df.diameter_m.min()

        # 4. open hole
        endCARFIN(LGR_NAME,
                    reopen_ID,
                    self.LGR_sizes_x, 
                    self.main_grd_min_k+1, 
                    self.min_grd_size,
                    self.no_of_layers_in_OB,
                    O)

        # 2. done
        O.close()

        # for qc
        print ('Output LGR CARFIN to: ', os.path.abspath(out_fname))

