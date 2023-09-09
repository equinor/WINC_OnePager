
import numpy as np
import pandas as pd

from src.WellClass.libs.utils.LGR_bbox import (
    get_ij_indices,
    get_k_indices,
)

class WellDataFrame:

    def __init__(self, my_well):

        # # Dataframes for drilling, casings, borehole and barriers

        # dataframes for the well
        self.drilling_df = pd.DataFrame(my_well.drilling)
        self.casings_df = pd.DataFrame(my_well.casings)
        self.borehole_df = pd.DataFrame(my_well.borehole)

        # and for the barriers
        self.barriers_df = pd.DataFrame(my_well.barriers)
        self.barriers_mod_df = pd.DataFrame(my_well.barriers_mod)


    def compute_bbox(self, mesh_df: pd.DataFrame, nxy: int):
        """ Compute bounding boxes for drillings, casings and barriers.
        """

        # for convenience
        casings_df = self.casings_df
        drilling_df = self.drilling_df
        barriers_df = self.barriers_df

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

        barriers_df['k_min'] = np.nan
        barriers_df['k_max'] = np.nan

        for idx, row in barriers_df.iterrows():
            
            top, bottom  = row['top_msl'], row['bottom_msl']

            # convert to indices
            k_min, k_max = get_k_indices(mesh_df, top, bottom)
            
            barriers_df.loc[idx, 'k_min'] = k_min
            barriers_df.loc[idx, 'k_max'] = k_max

    # def compute_annulus_thickness(self):
    #     """ Compute annulus thickness. For simplicty it assumes 
    #         annulus between casing and openhole as entire annulus
    #     """
    #     # for convenience
    #     casings_df = self.casings_df
    #     drilling_df = self.drilling_df

    #     #
    #     casings_df['ann_od_m'] = np.nan
    #     casings_df['ann_bottom_msl'] = np.nan

    #     # casing
    #     for idx, row in casings_df[::-1].iterrows():

    #         # get id from casing
    #         d, top, bottom = row[['diameter_m', 'top_msl', 'bottom_msl']]

    #         # get od from drilling
    #         hole = drilling_df[drilling_df['diameter_m'] > d].iloc[-1]
        
    #         hole_top, hole_bottom, hole_d = hole[['top_msl', 'bottom_msl', 'diameter_m']]

    #         casings_df.loc[idx, 'ann_od_m'] = hole_d       # outer diameter in meters
    #         casings_df.loc[idx, 'ann_bottom_msl'] = hole_bottom

    #     #Compute inner area
    #     casings_df['A_i'] = np.pi * (casings_df['diameter_m']/2)**2

    #     #Compute outer area
    #     casings_df['A_o'] = np.pi * (casings_df['ann_od_m']/2)**2

    #     # annulus thickness: (od-id)/2
    #     casings_df['thick_m'] = (casings_df['ann_od_m'] - casings_df['diameter_m'])/2
