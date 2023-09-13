
import numpy as np
import pandas as pd

def compute_annulus(casings: dict, drilling: dict) -> dict:
    """ Compute annulus fields. For simplicty it assumes 
        annulus between casing and openhole as entire annulus

        Args:
            casings (dict): contains casing information
            drilling (dict): contains drilling information

        Returns:
            annulus (dict): contains annulus information
    """

    casings_df  = pd.DataFrame(casings)
    drilling_df = pd.DataFrame(drilling)

    annulus_fields = []
    # casing
    for idx, row in casings_df[::-1].iterrows():

        # get id from casing
        d, top, bottom = row[['diameter_m', 'top_msl', 'bottom_msl']]

        # get od from drilling
        hole = drilling_df[drilling_df['diameter_m'] > d].iloc[-1]
    
        # extract the fields
        hole_top, hole_bottom, hole_d = hole[['top_msl', 'bottom_msl', 'diameter_m']]

        # collect them
        annulus_fields.append((hole_d, hole_top, hole_bottom))

    # make a dataframe
    annulus_df = pd.DataFrame(data=annulus_fields, columns=['ann_od_m', 'top_msl', 'bottom_msl'])

    #Compute inner area
    annulus_df['A_i'] = np.pi * (casings_df['diameter_m']/2)**2

    #Compute outer area
    annulus_df['A_o'] = np.pi * (annulus_df['ann_od_m']/2)**2

    # annulus thickness: (od-id)/2
    annulus_df['thick_m'] = (annulus_df['ann_od_m'] - casings_df['diameter_m'])/2

    return annulus_df.to_dict()

