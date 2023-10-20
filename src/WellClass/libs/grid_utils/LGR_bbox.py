

# handle type hints problem for python version < 3.10
from typing import Union, Tuple

import pandas as pd

def get_k_indices(df: pd.DataFrame, top: Union[float, int], bottom: Union[float, int]) -> Tuple[int, int]:
        """
        Takes the mesh data frame and a value of top and bottom depth interval.
        
        Args:
            df (pd.DataFrame): dataframe
            top (float): top depth value of the well element
            bottom (float): bottom depth value of the well element

        Returns:
            tuple: the min and max k indices of well element
        """
    
        # k_min
        if top <= df['Zcorn_top'].min():
            
                k_min = df['k'].min()
            
        elif top in df['Zcorn_top'].values:
            
                k_min = df.query('Zcorn_top==@top & Zcorn_bottom>=@top')['k'].iloc[0]
        else:
                k_min = df.query('Zcorn_top<=@top & Zcorn_bottom>=@top')['k'].iloc[0]

        # k_max
        if bottom >= df['Zcorn_bottom'].max():
            
                k_max = df['k'].max()

        elif bottom in df['Zcorn_bottom'].values:
                # print(bottom)
                k_max = df.query('Zcorn_top<=@bottom & Zcorn_bottom==@bottom')['k'].iloc[0]
        else:
                k_max = df.query('Zcorn_top<=@bottom & Zcorn_bottom>=@bottom')['k'].iloc[0]


        return k_min, k_max

def get_ij_indices(nxy: int, n_grd: int) -> Tuple[int, int]:
    """ compute x-y min/max indices

        Args:
            nxy (int): total grid size of x-y finer grid, i.e., refined x-y size of the center coarse grid
            n_grd (int): number of x grid of that well element, i.e., row['n_grd_id']

        Returns:
            tuple: min/max indices of given well element
    """
    # x-y ranges
    ij_min = (nxy - n_grd)//2
    ij_max = ij_min + n_grd - 1
    
    return ij_min, ij_max

def compute_bbox(mesh_df: pd.DataFrame, 
                 section_df: pd.DataFrame,
                 *, 
                 nxy: int,
                 is_casing: bool=False,
                 maxDepth: float=None) -> None:
        """ Compute bounding boxes for drillings, casings and barriers

            Args:

                mesh_df (pd.DataFrame): dataframe for refined grid
                section_df (pd.DataFrame): dataframe to be filled with bbox
                nxy (int): total grid size of x-y finer grid, i.e., refined x-y size of the center coarse grid
                is_casing (bool): boolean flag, only True for casing section
                maxDepth (float): maximum depth
        """

        # z
        section_df['k_min'] = -1
        section_df['k_max'] = -1
        # xy
        section_df['ij_min'] = -1
        section_df['ij_max'] = -1

        # cement bond, only for casings
        if is_casing:
                # z
                section_df['toc_k_min'] = -1
                section_df['toc_k_max'] = -1

        for idx, row in section_df.iterrows():
                
                # A) z ranges
                top, bottom = row['top_msl'], row['bottom_msl']

                if maxDepth and top > maxDepth:
                        continue
                
                # k ranges
                k_min, k_max = get_k_indices(mesh_df, top, bottom)

                # B) xy ranges
                ij_min, ij_max = get_ij_indices(nxy, row['n_grd_id'])

                # to dataframe
                section_df.loc[idx, 'k_min'] = k_min
                section_df.loc[idx, 'k_max'] = k_max
                section_df.loc[idx, 'ij_min'] = ij_min
                section_df.loc[idx, 'ij_max'] = ij_max

                # C) cement bond
                if is_casing:

                        # z ranges
                        toc, boc = row['toc_msl'], row['boc_msl']

                        # convert to indices
                        toc_k_min, toc_k_max = get_k_indices(mesh_df, toc, boc)

                        # to dataframe
                        section_df.loc[idx, 'toc_k_min'] = toc_k_min
                        section_df.loc[idx, 'toc_k_max'] = toc_k_max
