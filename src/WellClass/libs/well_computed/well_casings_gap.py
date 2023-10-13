import pandas as pd

def gap_casings(drilling_df: pd.DataFrame, 
                casings_df: pd.DataFrame):
    """ casings for GaP
    """

    # drilling columns
    drilling_columns = ['diameter_m', 'top_msl', 'bottom_msl']
    # casing columns
    casing_columns = ['diameter_m', 'top_msl', 'bottom_msl', 'toc_msl', 'boc_msl', 'cb_perm']

    # build dataframe
    new_casing_list = []
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
                # cement bond
                toc_msl = max(row_casing_df['toc_msl'], top_msl)
                boc_msl = min(row_casing_df['boc_msl'], bot_msl)

                # save them to list
                new_casing_list.append((row_casing_df['diameter_m'],
                                        top_msl,
                                        bot_msl,
                                        toc_msl,
                                        boc_msl,
                                        row_casing_df['cb_perm']))

        except Exception as e:

            break

    # build dataframe
    new_casing_df = pd.DataFrame(data=new_casing_list, 
                                 columns=casing_columns)

    return new_casing_df
