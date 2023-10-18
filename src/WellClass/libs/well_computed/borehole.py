
import pandas as pd

def compute_borehole(casings: dict, drilling: dict) -> dict:
    '''
    Routine to compute the effective open borehole. 
    Borehole relies mainly on the casings. It comes to drilling only for the places where casings are not available

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

    # saved for next section
    top_msl_saved = -1.0
    diameter_m_saved = -1.0

    # collect open holes
    borehole_list = []
    for ic, idx in enumerate(casings_df.index):

        # extract casings data
        row_casing_df = casings_df.loc[idx, casing_columns]

        # collect fields
            
        # casing info
        if ic == 0:
            top_msl = row_casing_df['top_msl']
        else:
            top_msl = top_msl_saved
        bot_msl = row_casing_df['bottom_msl']
        diameter_m = row_casing_df['diameter_m']

        # save them to list
        borehole_list.append((top_msl,
                              bot_msl, 
                              diameter_m
                                ))

        # for next casing
        top_msl_saved = bot_msl
        diameter_m_saved = diameter_m

    # TODO(hzh): a trade off, extend open hole in the last casing to bottom of corresponding drilling section
    # handle the transistion zone
    for idx in drilling_df.index:

        # extract drilling data
        row_drilling_df = drilling_df.loc[idx, drilling_columns]

        # locate the last casing
        if row_drilling_df['bottom_msl'] >= top_msl_saved:

            # the last section
            bot_msl = row_drilling_df['bottom_msl']

            # locate the last section
            last_one = borehole_list[-1]
            assert last_one[1] == top_msl_saved and last_one[-1] == diameter_m_saved

            # replace it
            borehole_list[-1] = (last_one[0], bot_msl, last_one[-1])

            # for next one
            top_msl_saved = bot_msl

            break

    # append the rest with drilling
    for ic, idx in enumerate(drilling_df.index):

        # extract drilling data
        row_drilling_df = drilling_df.loc[idx, drilling_columns]

        if row_drilling_df['bottom_msl'] > top_msl_saved:
            # save them to list
            borehole_list.append((row_drilling_df['top_msl'],
                                    row_drilling_df['bottom_msl'], 
                                    row_drilling_df['diameter_m'] 
                                    ))
    # build dataframe
    borehole_df = pd.DataFrame(data=borehole_list, columns=['top_msl', 'bottom_msl', 'id_m'])

    return borehole_df.to_dict()

