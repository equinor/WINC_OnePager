
import pandas as pd

def get_k_indices(df: pd.DataFrame, top: float|int, bottom: float|int) -> tuple[int, int]:
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

def get_ij_indices(nx: int, n_grd: int) -> tuple[int, int]:
    """ compute x-y min/max indices

        Args:
            nx (int): x grid size of finer grid
            n_grd (int): number of x grid of that well element, i.e., row['n_grd_id']

        Returns:
            tuple: min/max indices of given well element
    """
    # x-y ranges
    ij_min = (nx - n_grd)//2
    ij_max = ij_min + n_grd - 1
    
    return ij_min, ij_max

