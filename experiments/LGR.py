

from src.WellClass.libs.utils.LGR_grid_utils import (
    compute_ngrd,
    generate_LGR_xy,
    generate_LGR_z,
)

from .grid_coarse import GridCoarse

class LGR:

    def __init__(self,
                 grid_init: GridCoarse, 
                 drilling_df, casings_df, 
                 Ali_way):
        """ LGR grid information in x, y, z directions 
        """
        # initialize grid parameters
        self.main_grd_dx, self.main_grd_dy = grid_init.main_grd_dx, grid_init.main_grd_dy

        # DZs for reservoir and overburden
        self.DZ_rsrv = grid_init.DZ_rsrv
        self.DZ_ovb_coarse = grid_init.DZ_ovb_coarse

        # reference depth wheer LGR starts
        ref_depth = 0
        if Ali_way: 
            ref_depth = grid_init.ref_depth

        # ### 1. Compute minimum grid size

        min_grd_size = self._compute_min_grd_size(self, casings_df)
        print(f'Minimimum grid size is {min_grd_size*100:.2f} cm')

        # TODO(hzh): manually set it
        if Ali_way:
            min_grd_size = 0.05

        self.min_grd_size = min_grd_size

        # 1. compute number of LGR grids for drilling, casing and borehole, respectively

        # for drilling
        # drilling_df['n_grd_id']  = drilling_df['diameter_m'].map(lambda x: compute_ngrd(x, min_grd_size))
        drilling_series = drilling_df['diameter_m'].map(lambda x: compute_ngrd(x, min_grd_size))

        # 2.1 Number of cells representing horizontal LGR
        self.no_latral_fine_grd = drilling_series.max()

        # ### 2. Compute LGR grid sizes in x-y directions


    def _compute_min_grd_size(self, casings_df):

        # 0. minimum grid size

        # minimum grid size depends on minimum annulus thickness
        min_grd_size = casings_df['thick_m'].min()

        if min_grd_size < 0.05:
            min_grd_size = 0.05

        return min_grd_size
    
    def _compute_LGR_sizes_xy(self, 
                              Ali_way):
        """ Compute LGR grid sizes in x-y directions
        """

        # 2.2 generate the LGR grid sizes in x-y
        LGR_sizes_x, LGR_sizes_y, _ = generate_LGR_xy(no_latral_fine_grd, 
                                                        min_grd_size, 
                                                        main_grd_dx, main_grd_dy,
                                                        Ali_way=Ali_way)

        return LGR_sizes_x, LGR_sizes_y
    
    def _compute_LGR_sizes_z(self, 
                             DZ_rsrv, DZ_ovb_coarse, 
                             ref_depth):
        """ Compute LGR grid sizes in z direction
        """
        # LGR_sizes_z, LGR_numb_z, LGR_depths, _ = generate_LGR_z(DZ_rsrv, DZ_ovb_coarse)
        # TODO(hzh): to make LGR starts at ref_depth
        LGR_sizes_z, LGR_numb_z, LGR_depths, _ = generate_LGR_z(DZ_rsrv, DZ_ovb_coarse, ref_depth)

        return LGR_sizes_z, LGR_numb_z, LGR_depths

