
import pandas as pd


def compute_borehole(casings: dict, drilling: dict) -> dict:
    '''
    Update in Routine to compute the effective open borehole. 
    The previous version remains stored as borehole_alt
    Borehole relies mainly on the casings. It comes to drilling only for the places where casings are not available

        Args:
            casings (dict): contains casing information
            drilling (dict): contains drilling information

        Returns:
            borehole (dict): contains borehole information
    '''
    casings_df  = pd.DataFrame(casings)
    drilling_df = pd.DataFrame(drilling)    
    
    #Merge casings and drilling tables
    well_concat = pd.concat([casings_df, drilling_df])

    #Sort new table, firstly by top_rkb, and  secondly by diameter
    well_concat = well_concat.sort_values(by=['top_msl', 'diameter_m'])
    well_concat.reset_index(inplace=True)

    #iterate over merged table
    borehole_list = []

    z_values = well_concat['top_msl'].drop_duplicates()

    #very large value to initiate iteration
    previous_diam = 999

    #iterate over top_msl values
    for z in z_values:
            #query rows in merged table that have the z value as top_msl
            query_df = well_concat.query('top_msl==@z')

            #retrieve minimum diameter
            min_diam_m = query_df['diameter_m'].min()

            #check if diameter is smaller tahn previous_diam
            if min_diam_m < previous_diam:
                    #append top_msl and diameter
                    borehole_list.append((z, min_diam_m))
                    
                    #update previous_diam for next iteration
                    previous_diam = min_diam_m


    #create dataframe
    borehole_df = pd.DataFrame(data=borehole_list, columns=['top_msl', 'diameter_m'])

    #shit top_msl to create bottom_msl column
    borehole_df['bottom_msl'] = borehole_df['top_msl'].shift(-1)

    # Add TD in MSL to last record of bottom_msl
    borehole_df.loc[borehole_df.index[-1], 'bottom_msl'] = well_concat['bottom_msl'].max()

    borehole_df.to_dict()


    return borehole_df.to_dict()

