
from typing import Tuple

import numpy as np
import pandas as pd

from .grid_coarse import GridCoarse

from .LGR_bbox import (
    get_ij_indices,
    get_k_indices,
)

class GridRefineBase:

    def __init__(self, 
                 grid_coarse: GridCoarse,
                 LGR_sizes_x, LGR_sizes_y, LGR_sizes_z,
                 ):
        """ dataframe for LGR mesh for the center coarse grid
        """

        # dx/dy of coarse grid
        self.main_grd_dx = grid_coarse.main_grd_dx
        self.main_grd_dy = grid_coarse.main_grd_dy

        # LGR dimensions
        nx = len(LGR_sizes_x)
        ny = len(LGR_sizes_y)
        nz = len(LGR_sizes_z)

        print(f'nx={nx}, ny={ny}, nz={nz}')

        # set dimensions
        self.nx, self.ny, self.nz = nx, ny, nz

        #Create i, j, k indices
        cell_ijk = np.indices((nx, ny, nz))
        cell_ijk = cell_ijk.reshape((3, nx * ny * nz)).T

        #Create LGR Dataframe with indices
        mesh_df = pd.DataFrame(data = cell_ijk, columns = ['i', 'j', 'k'])
        #
        mesh_df.sort_values(by=['k', 'i'], inplace = True)
        mesh_df.reset_index(inplace=True, drop=True)

        #
        self.mesh_df = mesh_df


        # set up fields
        self._set_cell_intervals(self.mesh_df, 
                                 LGR_sizes_x, 
                                 LGR_sizes_z, 
                                 LGR_sizes_y)
        self._set_cell_coords(self.mesh_df)
        self._set_z_corners(self.mesh_df)       

        # TODO(hzh): do we suppose to use mid index of coarse grid?

        # indices for center finer grid
        mid_i = grid_coarse.main_grd_i       # mesh_df.i.max()//2
        mid_j = grid_coarse.main_grd_j       # mesh_df.j.max()//2

        self._upscale_properties(grid_coarse.grid_init, mid_i, mid_j)
        
    def _set_cell_intervals(self, mesh_df, 
                            LGR_sizes_x, 
                            LGR_sizes_z, 
                            LGR_sizes_y):
        """ Creat DX, DY, DZ for LGR mesh
        """
 
        # mesh
        DX_grid, DZ_grid, DY_grid = np.meshgrid(LGR_sizes_x, LGR_sizes_z, LGR_sizes_y)

        # internals
        mesh_df['DX'] = DX_grid.flatten()
        mesh_df['DY'] = DY_grid.flatten()
        mesh_df['DZ'] = DZ_grid.flatten()

    def _set_cell_coords(self, mesh_df):
        """ Create cell coordinate X, Y, Z for LGR mesh
        """

        # cell coordinates
        xcoord = (mesh_df.query("j==0&k==0").DX.cumsum() - mesh_df.query("j==0&k==0").DX/2).values
        ycoord = (mesh_df.query("i==0&k==0").DY.cumsum() - mesh_df.query("i==0&k==0").DY/2).values
        zcoord = (mesh_df.query("i==0&j==0").DZ.cumsum() - mesh_df.query("i==0&j==0").DZ/2).values

        # TODO(hzh): a bug?
        # map_X = dict(zip(mesh_df.query("j==0&j==0")['i'], xcoord))
        map_X = dict(zip(mesh_df.query("j==0&k==0")['i'], xcoord))
        map_Y = dict(zip(mesh_df.query("i==0&k==0")['j'], ycoord))
        map_Z = dict(zip(mesh_df.query("i==0&j==0")['k'], zcoord))

        # save cell coordinates X, Y, Z to dataframe
        mesh_df['X'] = mesh_df['i'].map(map_X)
        mesh_df['Y'] = mesh_df['j'].map(map_Y)
        mesh_df['Z'] = mesh_df['k'].map(map_Z)

    def _set_z_corners(self, mesh_df):
        """ save Corner Z points to dataframe 
        """

        mesh_df['Zcorn_top'] = mesh_df['Z'] - mesh_df['DZ']/2
        mesh_df['Zcorn_bottom'] = mesh_df['Z'] + mesh_df['DZ']/2

    def _upscale_properties(self, grid_init, mid_i, mid_j):
        """ Upscale coarse properties to LGR mesh
        """

        # for convenience only
        mesh_df = self.mesh_df

        # properties
        fields = ['PORV', 
                  'PERMX', 'PERMY', 'PERMZ', 
                  'MULTX', 'MULTY', 'MULTZ', 
                  'MULTX-', 'MULTY-', 'MULTZ-', 
                  'PORO']

        # Upscale coarse properties to LGR grids

        for field in fields:
            # mesh_df[field] = np.nan
            mesh_df[field] = 0.0                    # TODO(hzh): what should I put here?
            
        for idx, row in grid_init.query('i==@mid_i & j==@mid_j').iterrows():

            # switch to corner coords, coarse grid
            top  = row.Z - row.DZ/2
            base = row.Z + row.DZ/2

            for field in fields:

                mesh_df.loc[(mesh_df['Z']>=top) & (mesh_df['Z']<base), field] = row[field]



    def _set_material_type(self,
                           drilling_df, 
                           casings_df, 
                           barriers_mod_df):
        """ Assign material types, such as openholes, overburden, cement bond, etc.
        """

        # only for convenience
        mesh_df = self.mesh_df

        # set default material to 'overburden'
        mesh_df['material'] = 'overburden'

        # ### 1. Drilling
        for idx, row in drilling_df.iterrows():
            
            top, bottom = row['top_msl'], row['bottom_msl']

            if top < mesh_df['Zcorn_bottom'].max():

                # extract bounding box
                k_min, k_max = row['k_min'], row['k_max']
                ij_min, ij_max = row['ij_min'], row['ij_max']
                
                # 1.1 set material type to openhole
                criteria =  '(k >= @k_min)  & (k <= @k_max) & \
                             (i >= @ij_min) & (i <= @ij_max) & \
                             (j >= @ij_min) & (j <= @ij_max)'
                mesh_df.loc[mesh_df.eval(criteria), 'material'] = 'openhole'

        # ### 2. Casings
        for ic, (idx, row) in enumerate(casings_df.iterrows()):

            # extract bounding box
            k_min, k_max = row['k_min'], row['k_max']
            ij_min, ij_max = row['ij_min'], row['ij_max']
            toc_k_min, toc_k_max = row['toc_k_min'], row['toc_k_max']
            
            # 2.1 set material type to annulus
            # x
            criteria_i =  '(material == "openhole") & \
                           (k >= @k_min) & (k <= @k_max) & \
                            ((i < @ij_min) | (i > @ij_max))'
            mesh_df.loc[mesh_df.eval(criteria_i), 'material'] = 'annulus'
            # y
            criteria_j =  '(material == "openhole") & \
                           (k >= @k_min) & (k <= @k_max) & \
                           ((j < @ij_min) | (j > @ij_max))'
            mesh_df.loc[mesh_df.eval(criteria_j), 'material'] = 'annulus'
            
            # 2.2 set material type to cement_bond
            criteria = '(material == "annulus") & \
                        (k >= @toc_k_min) & (k <= @toc_k_max)' 
            mesh_df.loc[mesh_df.eval(criteria), 'material'] = f'cement_bond_{ic}'

            # 2.3 set material type to openhole
            criteria = '(material == "annulus")'  
            mesh_df.loc[mesh_df.eval(criteria), 'material'] = 'openhole'
            # mesh_df.loc[mesh_df.eval(criteria_j), 'material'] = 'cementbond'   

        # ### 3. Barriers
        for ib, (idx, row) in enumerate(barriers_mod_df.iterrows()):
            
            b_k_min, b_k_max = row['k_min'], row['k_max']
            
            criteria = '(material == "openhole") & \
                        (k >= @b_k_min) & (k <= @b_k_max)' 
            mesh_df.loc[mesh_df.eval(criteria), 'material'] = f'barrier_{ib}'


    def _set_permeability(self, drilling_df, casings_df, barriers_mod_df):
        """ Actual function to assign permeability according to material type
        """

        # for convenience only
        mesh_df = self.mesh_df

        # set permeability according to material type

        # 1. openhole
        oh_perm = drilling_df['oh_perm'].iloc[0]
        criteria = 'material == "openhole"'
        mesh_df.loc[mesh_df.eval(criteria), 'PERMX'] = oh_perm

        # 2. cement bond
        for ic, (_, row) in enumerate(casings_df.iterrows()):
            cb_perm = row['cb_perm']
            criteria = f'material == "cement_bond_{ic}"'
            mesh_df.loc[mesh_df.eval(criteria), 'PERMX'] = cb_perm

        # 3. barrier
        for ib, (_, row) in enumerate(barriers_mod_df.iterrows()):
            barrier_perm = row['barrier_perm']
            criteria = f'material == "barrier_{ib}"'
            mesh_df.loc[mesh_df.eval(criteria), 'PERMX'] = barrier_perm

    def _compute_bbox(self, drilling_df, casings_df, barriers_mod_df):
        """ Compute bounding boxes for drillings, casings and barriers.
        """

        # for convenience
        mesh_df = self.mesh_df
        nxy = self.nx

        # ### 1. Drillings

        drilling_df['k_min'] = np.nan
        drilling_df['k_max'] = np.nan
        drilling_df['ij_min'] = np.nan
        drilling_df['ij_max'] = np.nan

        for idx, row in drilling_df.iterrows():
            
            top, bottom = row['top_msl'], row['bottom_msl']

            if top < mesh_df['Zcorn_bottom'].max():
                
                # k ranges
                k_min, k_max = get_k_indices(mesh_df, top, bottom)

                # x-y ranges
                ij_min, ij_max = get_ij_indices(nxy, row['n_grd_id'])

                # to dataframe
                drilling_df.loc[idx, 'k_min'] = k_min
                drilling_df.loc[idx, 'k_max'] = k_max
                drilling_df.loc[idx, 'ij_min'] = ij_min
                drilling_df.loc[idx, 'ij_max'] = ij_max# # Bounding box for well elements

        # ### 2. Casings

        # casing, k
        casings_df['k_min'] = np.nan
        casings_df['k_max'] = np.nan
        # cement bond, k
        casings_df['toc_k_min'] = np.nan
        casings_df['toc_k_max'] = np.nan
        # casing, xy
        casings_df['ij_min'] = np.nan
        casings_df['ij_max'] = np.nan

        for idx, row in casings_df.iterrows():

            # A) casing, z ranges
            top, bottom  = row['top_msl'], row['bottom_msl']
            
            # convert to indices
            k_min, k_max = get_k_indices(mesh_df, top, bottom)
            
            # B) cement, z ranges
            toc, boc = row['toc_msl'], row['boc_msl']

            # convert to indices
            toc_k_min, toc_k_max = get_k_indices(mesh_df, toc, boc)
            
            # C) xy ranges
            ij_min, ij_max = get_ij_indices(nxy, row['n_grd_id'])
            
            # to dataframe
            casings_df.loc[idx, 'k_min'] = k_min
            casings_df.loc[idx, 'k_max'] = k_max
            casings_df.loc[idx, 'toc_k_min'] = toc_k_min
            casings_df.loc[idx, 'toc_k_max'] = toc_k_max
            casings_df.loc[idx, 'ij_min'] = ij_min
            casings_df.loc[idx, 'ij_max'] = ij_max

        # ### 3. Barriers

        barriers_mod_df['k_min'] = np.nan
        barriers_mod_df['k_max'] = np.nan

        for idx, row in barriers_mod_df.iterrows():
            
            top, bottom  = row['top_msl'], row['bottom_msl']

            # convert to indices
            k_min, k_max = get_k_indices(mesh_df, top, bottom)
            
            barriers_mod_df.loc[idx, 'k_min'] = k_min
            barriers_mod_df.loc[idx, 'k_max'] = k_max