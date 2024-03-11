
from collections import defaultdict

import pandas as pd

def compute_barriers_diam(barriers: dict, borehole: dict) -> dict:
    '''
    Processes barriers data. Takes the cumputed borehole data to recompute 
    the barrier intervals in cases where the barriers are sitting along differen hole sizes
    '''

    barriers_df = pd.DataFrame(barriers)
    borehole_df = pd.DataFrame(borehole)

    barriers_fmt = defaultdict(lambda: defaultdict())

    for idx, row in barriers_df.iterrows():

        barr_idx = idx

        barrier = row['barrier_name']   # row[0]
        # TODO(hzh): not used
        # b_type = row['barrier_type']    # row[1]
        top = row['top_msl']
        bottom = row['bottom_msl']

        # TODO(hzh): for compatibility with .csv file
        try:
            barrier_perm = row['barrier_perm']           # barrier permeability
        except KeyError:
            barrier_perm = 0.5    # make a fake one if doesn't exist

        check_top_int_top = top >= borehole_df['top_msl']
        check_bottom_int_top = top <= borehole_df['bottom_msl']

        check_top_int_bottom = bottom >= borehole_df['top_msl']
        check_bottom_int_bottom = bottom <= borehole_df['bottom_msl']

        barrier_topD = borehole_df[check_top_int_top & check_bottom_int_top]['id_m'].iloc[0]
        barrier_bottomD = borehole_df[check_top_int_bottom & check_bottom_int_bottom]['id_m'].iloc[0]

        if barrier_topD == barrier_bottomD:
            barrier_name = '{:s}_{:d}'.format(barrier, barr_idx)
            barriers_fmt[barrier_name]['b_name'] = barrier
            barriers_fmt[barrier_name]['top_msl'] = top
            barriers_fmt[barrier_name]['bottom_msl'] = bottom
            barriers_fmt[barrier_name]['diameter_m'] = barrier_topD
            barriers_fmt[barrier_name]['barrier_perm'] = barrier_perm
        else:
            barrier_query = borehole_df.query('id_m <= @barrier_topD and id_m >= @barrier_bottomD')
            for idx, sub_section in barrier_query.iterrows():
                barrier_name = '{:s}_{:d}'.format(barrier, barr_idx)
                barriers_fmt[barrier_name]['b_name'] = barrier

                barriers_fmt[barrier_name]['top_msl'] = max(sub_section['top_msl'], top)
                barriers_fmt[barrier_name]['bottom_msl'] = min(sub_section['bottom_msl'], bottom)
                barriers_fmt[barrier_name]['diameter_m'] = sub_section['id_m']
                barriers_fmt[barrier_name]['barrier_perm'] = barrier_perm
                barr_idx += 1

    barriers_fmt_df = pd.DataFrame(barriers_fmt).T
    
    return barriers_fmt_df.to_dict()


