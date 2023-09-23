
from typing import Tuple

import numpy as np
import pandas as pd

def extract_xz_corn_coords(grid_df: pd.DataFrame, 
                           sDX: float, 
                           sDY: float=None) -> Tuple[np.ndarray, np.ndarray]:
    """ extract grid x-z coordinates and shift them to the middle

        Args:
            grid_df (pd.DataFrame): grid dataframe
            sDX (float): the x distance to shift
            sDY (float): the y distance to shift

        Returns:

            Tuple[np.ndarray, np.ndarray]: The x and z coordinates

    """
    xcorn  = (grid_df.query("j==0&k==0").DX.cumsum()).values
    ycorn  = (grid_df.query("i==0&k==0").DY.cumsum()).values
    zcorn  = (grid_df.query("i==0&j==0").DZ.cumsum()).values

    # add origin coordinates
    xcorn = np.append(0, xcorn)
    ycorn = np.append(0, ycorn)
    zcorn = np.append(0, zcorn)

    # shift grid coordinates half-length in x-y directions
    # but not in z direction
    xcorn -= sDX
    if sDY:
        ycorn -= sDY

    return xcorn, zcorn

def extract_xz_prop_slice(grid_df: pd.DataFrame,
                          *, 
                          mid_j: int=None,
                          prop='PERMX') -> np.ndarray:
    """ generate x-z PERM slice

        Args:
            grid_df (pd.DataFrame): grid dataframe
            mid_j (int): given j index to extract slice
            prop (str): the proper name, default: PERMX

        Returns:
            np.ndarray: the property x-z slide
    """

    # the dimensions of the grid
    nx = grid_df.i.max() + 1
    ny = grid_df.j.max() + 1
    kmin = grid_df.k.min()
    kmax = grid_df.k.max()
    nz = kmax - kmin + 1

    # extract only mid j
    if mid_j is None:
        mid_j = ny//2

    # extract 2D xz slice at middle of y
    XZ_slice = grid_df.query('j==@mid_j')

    # extract permeability
    Z = XZ_slice[prop].values.reshape(nz, nx)

    return Z
