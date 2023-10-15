
import numpy as np
import pandas as pd

def compute_borehole(casings: dict, drilling: dict) -> dict:
    '''
    Routine to compute the effective open borehole. 
    It should be the casing, or the drilling if no casing is available

        Args:
            casings (dict): contains casing information
            drilling (dict): contains drilling information

        Returns:
            borehole (dict): contains borehole information
    '''
    
    casings_df  = pd.DataFrame(casings)
    drilling_df = pd.DataFrame(drilling)

    # drilling columns
    drilling_columns = ['diameter_m', 'top_msl', 'bottom_msl']
    # casing columns
    casing_columns = ['diameter_m', 'top_msl', 'bottom_msl']

    # build dataframe
    borehole_list = []
    for idx in drilling_df.index:

        # extract drilling data
        row_drilling_df = drilling_df.loc[idx, drilling_columns]

        try:
            # extract casings data
            row_casing_df = casings_df.loc[idx, casing_columns]

            # collect fields
            if row_casing_df['bottom_msl'] < row_drilling_df['bottom_msl']:
                
                # casing
                top_msl = max(row_drilling_df['top_msl'], row_casing_df['top_msl'])
                bot_msl = min(row_drilling_df['bottom_msl'], row_casing_df['bottom_msl'])

                # save them to list
                borehole_list.append((top_msl,
                                        bot_msl, 
                                        row_casing_df['diameter_m']
                                        ))

        except Exception as e:

                # save them to list
                borehole_list.append((row_drilling_df['top_msl'],
                                      row_drilling_df['bottom_msl'], 
                                      row_drilling_df['diameter_m'] 
                                      ))


    borehole_df = pd.DataFrame(data=borehole_list, columns=['top_msl', 'bottom_msl', 'id_m'])

    return borehole_df.to_dict()

