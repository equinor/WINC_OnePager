
from src.WellClass.libs.utils.LGR_grid_utils import (
    compute_ngrd,
    generate_LGR_xy,
    generate_LGR_z,
)

class LGR:

    def __init__(self, drilling_df, casings_df, Ali_way):
        """ LGR grid information in x, y, z directions 
        """

        # ### 1. Compute minimum grid size

        # 0. minimum grid size

        # minimum grid size depends on minimum annulus thickness
        min_grd_size = casings_df['thick_m'].min()

        if min_grd_size < 0.05:
            min_grd_size = 0.05

        print(f'Minimimum grid size is {min_grd_size*100:.2f} cm')

        # TODO(hzh): manually set it
        if Ali_way:
            min_grd_size = 0.05

        # 1. compute number of LGR grids for drilling, casing and borehole, respectively

        # for drilling
        # drilling_df['n_grd_id']  = drilling_df['diameter_m'].map(lambda x: compute_ngrd(x, min_grd_size))
        drilling_series = drilling_df['diameter_m'].map(lambda x: compute_ngrd(x, min_grd_size))

        # ### 2. Compute LGR grid sizes in x-y directions

        # 2.1 Number of cells representing horizontal LGR
        no_latral_fine_grd = drilling_series.max()

        # 2.2 generate the LGR grid sizes in x-y
        LGR_sizes_x, LGR_sizes_y, _ = generate_LGR_xy(no_latral_fine_grd, 
                                                        min_grd_size, 
                                                        main_grd_dx, main_grd_dy,
                                                        Ali_way=Ali_way)

        # ### 3. Compute LGR grid sizes in z direction

        # LGR_sizes_z, LGR_numb_z, LGR_depths, _ = generate_LGR_z(DZ_rsrv, DZ_ovb_coarse)
        # TODO(hzh): to make LGR starts at ref_depth
        LGR_sizes_z, LGR_numb_z, LGR_depths, _ = generate_LGR_z(DZ_rsrv, DZ_ovb_coarse, ref_depth)
