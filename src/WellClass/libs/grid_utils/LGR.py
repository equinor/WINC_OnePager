
import pandas as pd

from .LGR_grid_utils import (
    compute_ngrd,
    generate_LGR_xy,
    generate_LGR_z,
)

from .grid_coarse import GridCoarse

from .df2gap import (
    to_gap_casing_list,
    to_gap_barrier_list
)

from src.GaP.libs.carfin import build_grdecl

class LGR:

    def __init__(self,
                 grid_init: GridCoarse, 
                 annulus_df: pd.DataFrame,
                 drilling_df: pd.DataFrame, 
                 casings_df: pd.DataFrame, 
                 borehole_df: pd.DataFrame,
                 Ali_way: bool):
        """ LGR grid information in x, y, z directions. We are going to compute the grid sizes in lateral (x and y) and vertical directions

            Args:
                grid_init (GridCoarse): all information about coarse grid
                annulus_df (pd.DataFrame): information about annulus
                drilling_df (pd.DataFrame): information about drilling
                casings_df (pd.DataFrame): information about casings and cement-bond
                borehold_df (pd.DataFrame): borehole
                Ali_way (bool): use Ali's algorithm to compute lateral grids and apply refdepth in z direction
        """

        # initialize coarse grid parameters
        self.NX, self.NY = grid_init.NX, grid_init.NY
        self.main_grd_dx, self.main_grd_dy = grid_init.main_grd_dx, grid_init.main_grd_dy
        self.main_grd_i, self.main_grd_j = grid_init.main_grd_i, grid_init.main_grd_j
        self.main_grd_min_k, self.main_grd_max_k = grid_init.main_grd_min_k, grid_init.main_grd_max_k

        # DZs for reservoir and overburden
        self.DZ_rsrv = grid_init.DZ_rsrv
        self.DZ_ovb_coarse = grid_init.DZ_ovb_coarse

        # number of layers of ovb
        self.no_of_layers_in_OB = grid_init.no_of_layers_in_OB

        # reference depth wheer LGR starts
        self.ref_depth = 0
        if Ali_way: 
            self.ref_depth = grid_init.ref_depth

        ######################################

        # ### 1. Compute minimum grid size

        min_grd_size = self._compute_min_grd_size(annulus_df, Ali_way)
        self.min_grd_size = min_grd_size

        # #### 2. Compute number of cells of horizontal LGR
        self.no_latral_fine_grd = self._compute_no_latral_fine_grd(drilling_df, casings_df, borehole_df)

        # #### 3. compute LGR sizes

        # 3.1 compute LGR sizes in x-y directions
        self.LGR_sizes_x, self.LGR_sizes_y = self._compute_LGR_sizes_xy(Ali_way)

        # 3.2 comptue LGR sizes in z direction
        self.LGR_sizes_z, self.LGR_numb_z, self.LGR_depths = self._compute_LGR_sizes_z()

    def _compute_min_grd_size(self, annulus_df, Ali_way):

        # 0. minimum grid size

        # minimum grid size depends on minimum annulus thickness
        min_grd_size = annulus_df['thick_m'].min()

        if min_grd_size < 0.05:
            min_grd_size = 0.05

        print(f'Minimimum grid size is {min_grd_size*100:.2f} cm')

        # TODO(hzh): manually set it
        if Ali_way:
            min_grd_size = 0.05

        return min_grd_size
    
    def _compute_no_latral_fine_grd(self, drilling_df, casings_df, borehole_df):
        """ compute number of LGR laternal grids
        """

        # only for convenience
        min_grd_size = self.min_grd_size

        # n_grd_id for well elements
        drilling_df['n_grd_id']  = drilling_df['diameter_m'].map(lambda x: compute_ngrd(x, min_grd_size))
        casings_df[ 'n_grd_id']  = casings_df['diameter_m'].map(lambda x: compute_ngrd(x, min_grd_size))
        borehole_df['n_grd_id'] = borehole_df['id_m'].map(lambda x: compute_ngrd(x, min_grd_size))

        # 
        drilling_series = drilling_df['diameter_m'].map(lambda x: compute_ngrd(x, min_grd_size))

        return drilling_series.max()
    
    def _compute_LGR_sizes_xy(self, Ali_way: bool):
        """ Compute LGR grid sizes in x-y directions
        """

        # for convenience
        no_latral_fine_grd = self.no_latral_fine_grd
        min_grd_size = self.min_grd_size
        main_grd_dx = self.main_grd_dx
        main_grd_dy = self.main_grd_dy

        # 3.1 generate the LGR grid sizes in x-y
        LGR_sizes_x, LGR_sizes_y, _ = generate_LGR_xy(no_latral_fine_grd, 
                                                        min_grd_size, 
                                                        main_grd_dx, main_grd_dy,
                                                        Ali_way=Ali_way)

        return LGR_sizes_x, LGR_sizes_y
    
    def _compute_LGR_sizes_z(self):
        """ Compute LGR grid sizes in z direction
        """
        # for convenience
        DZ_rsrv = self.DZ_rsrv
        DZ_ovb_coarse = self.DZ_ovb_coarse 
        ref_depth = self.ref_depth

        # 3.2 generate the LGR grid sizes in x-y

        # TODO(hzh): to make LGR starts at ref_depth
        LGR_sizes_z, LGR_numb_z, LGR_depths, _ = generate_LGR_z(DZ_rsrv, DZ_ovb_coarse, ref_depth)

        return LGR_sizes_z, LGR_numb_z, LGR_depths

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


