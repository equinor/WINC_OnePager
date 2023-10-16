
import pandas as pd

def trim_casings_cement(casings_df: pd.DataFrame) -> pd.DataFrame:
    '''
    Routine to compute non-overlapped casings and cement bond 

        Args:

            casings_df (pd.DataFrame): information about casings and cement-bond

        Returns:
            trimmed casing information ready for GaP code
    '''

    # casing columns
    casing_columns = ['diameter_m', 'top_msl', 'bottom_msl', 'toc_msl', 'boc_msl', 'cb_perm']

    # saved for next section
    top_msl_saved = -1.0

    # collect open holes
    casing_cement_list = []
    for ic, idx in enumerate(casings_df.index):

        # extract casings data
        row_casing_df = casings_df.loc[idx, casing_columns]

        # collect fields
            
        # casing
        if ic == 0:
            top_msl = row_casing_df['top_msl']
        else:
            top_msl = top_msl_saved
        bot_msl = row_casing_df['bottom_msl']
        diameter_m = row_casing_df['diameter_m']

        # cement bond
        toc_msl = max(row_casing_df['toc_msl'], top_msl)
        boc_msl = min(row_casing_df['boc_msl'], bot_msl)

        # save them to list
        casing_cement_list.append((diameter_m,   # casing
                                    top_msl,
                                    bot_msl, 
                                    toc_msl,     # cement-bond
                                    boc_msl,
                                    row_casing_df['cb_perm']))

        # for next casing
        top_msl_saved = bot_msl

    # build dataframe
    casing_cement_df = pd.DataFrame(data=casing_cement_list, columns=casing_columns)

    return casing_cement_df

