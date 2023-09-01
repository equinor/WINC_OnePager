
import pandas as pd

def compute_cement_bond(casings: dict, drilling: dict) -> dict:
    '''
    Processes cement bond intervals. Reads both casing and borehole sizes to estimate cement bond width
    '''

    casings_df = pd.DataFrame(casings)
    drilling_df = pd.DataFrame(drilling)

    cb_fields = []

    for idx, row in casings_df[::-1].iterrows():
        
        d, top, bottom = row[['diameter_m', 'toc_msl', 'boc_msl']]

        hole = drilling_df[drilling_df['diameter_m'] > d].iloc[-1]
        hole_top, hole_bottom, hole_d = hole[['top_msl', 'bottom_msl', 'diameter_m']]

        if top >= hole_top and bottom <= hole_bottom:
            cb_fields.append((top, bottom, d, hole_d))
        else:
            bond_query = casings_df.query('diameter_m > @d & bottom_msl > @top')
            bond_query = bond_query.sort_values(by='diameter_m')
            temp_top_cb = bond_query['bottom_msl'].max()
            cb_fields.append((temp_top_cb, bottom, d, hole_d))

            for idx, row in bond_query.iterrows():
                if row['top_msl'] <= top:
                    section_top = top
                    section_bottom = temp_top_cb
                    cb_fields.append((section_top, section_bottom, d, row['diameter_m']))
                    break
                else:
                    section_top = row['top_msl']
                    section_bottom = temp_top_cb
                    cb_fields.append((section_top, section_bottom, d, row['diameter_m']))
                    temp_top_cb = row['top_msl']

    cement_bond_df = pd.DataFrame(data=cb_fields, columns=['top_msl', 'bottom_msl', 'id_m', 'od_m'])

    return cement_bond_df.to_dict()

