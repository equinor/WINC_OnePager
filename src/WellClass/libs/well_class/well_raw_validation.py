
""" validate input drilling, casing, cement bond, and barriers
"""
import pandas as pd

def valid_drilling(drilling_df: pd.DataFrame):
    """ validate drilling input
    """
    
    # ensure it is monotonic decreasing
    assert drilling_df['diameter_m'].is_monotonic_decreasing, "drilling needs to be sorted"

    # extract data
    top_bot = drilling_df.loc[:, ['top_msl', 'bottom_msl']].values
    top = top_bot[:, 0]
    bot = top_bot[:, 1]

    # ensure the depth order is correct
    assert (top < bot).all()

    # ensure the bottom depth matches its next top depth
    assert (top[1:] == bot[:-1]).all()

def valid_casings(casings_df: pd.DataFrame):
    """ validate casings input
    """

    # ensure it is monotonic decreasing
    assert casings_df['diameter_m'].is_monotonic_decreasing, "casings needs to be sorted"

    # extract data
    columns = ['diameter_m', 'top_msl', 'bottom_msl', 'toc_msl', 'boc_msl']
    casing_cols = casings_df.loc[:, columns]

    # ensure the depth order is correct
    assert (casing_cols['top_msl'].values < casing_cols['bottom_msl'].values).all()
    assert (casing_cols['toc_msl'].values < casing_cols['boc_msl'].values).all()

    # ensure the cement bond is within its corresponding casing
    assert (casing_cols['top_msl'].values <= casing_cols['toc_msl'].values).all() 
    assert (casing_cols['bottom_msl'].values >= casing_cols['boc_msl'].values).all()

def valid_barriers(barriers_df: pd.DataFrame):
    """ validate barriers input
    """
    pass
