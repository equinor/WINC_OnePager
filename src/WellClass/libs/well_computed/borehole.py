
import numpy as np
import pandas as pd

def compute_borehole(casings: dict, drilling: dict) -> dict:
    '''
    Routine to compute the effective open borehole. 
    It takes the original hole and substracts the diameter taken by cement bond.

        Args:
            casings (dict): contains casing information
            drilling (dict): contains drilling information

        Returns:
            borehold (dict): contains borehole information
    '''
    
    casings_df  = pd.DataFrame(casings)
    drilling_df = pd.DataFrame(drilling)

    well_concat = pd.concat([casings_df, drilling_df])
    top_z = well_concat['top_msl'].drop_duplicates().dropna().values
    top_z.sort()

    diam = 1 + well_concat['diameter_m'].max()

    borehole = []
    for z_value in top_z:
        
        z_query = well_concat.query('top_msl == @z_value')
        min_q = z_query.query('diameter_m == diameter_m.min()')
        if min_q.iloc[0]['diameter_m'] < diam:
            borehole.append(min_q.iloc[0][['top_msl', 'bottom_msl', 'diameter_m']].values.tolist())
            diam = min_q.iloc[0]['diameter_m']

    borehole = np.array(borehole)
    borehole[:-1, 1] = borehole[1:, 0]
    # borehole[1:, 0] = borehole[:-1, 1]

    borehole_df = pd.DataFrame(data=borehole, columns=['top_msl', 'bottom_msl', 'id_m'])

    return borehole_df.to_dict()

