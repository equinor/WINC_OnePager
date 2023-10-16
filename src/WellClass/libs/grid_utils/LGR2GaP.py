
from typing import TextIO

import pandas as pd

from src.GaP.libs.carfin.CARFIN_pipe_with_oph import CARFIN_pipe_with_oph
from src.GaP.libs.carfin.CARFIN_oph import CARFIN_oph
from src.GaP.libs.carfin.CARFIN_cement_bond import CARFIN_cement_bond
from src.GaP.libs.carfin.CARFIN_barrier import CARFIN_barrier

def df_to_gap_casing(drilling_df: pd.DataFrame, 
                     casings_df: pd.DataFrame,
                     LGR_NAME: str,
                     O: TextIO) -> None:
    """ convert casing dataframe to gap format

        Args:
            drilling_df (pd.DataFrame): dataframe for drilling
            casings_df (pd.DataFrame): dataframe for casings
            LGR_NAME (str): LGR name
            O (TextIO): opened file handle
    """
    # assume oh_perm is the same for open hole
    oh_perm = drilling_df['oh_perm'].values[0]
    drilling_k_maxs = drilling_df['oh_perm'].values

    # saved for open hole section
    k_max_oph_saved = -1

    # 1. casings
    for ic, idx in enumerate(casings_df.index):

        # bbox for casings
        pipe_columns = ['diameter_m', 'ij_min', 'ij_max', 'k_min', 'k_max']
        ID_pipe, ij_min_pipe, ij_max_pipe, k_min_pipe, k_max_pipe  = casings_df.loc[idx, pipe_columns]

        k_min_oph = k_min_pipe
        k_max_oph = k_max_pipe
        if ic == len(casings_df)-1:  # the last row
            # locate the nearest drilling k_max that is larger than k_max_pipe
            for kmax in drilling_k_maxs:
                if kmax >= k_max_pipe:
                    k_max_oph = k_max_oph_saved = kmax
                    break

        # 1.1 CARFIN output
        CARFIN_pipe_with_oph(ID_pipe,
                                int(ij_min_pipe+1), int(ij_max_pipe+1),
                                int(ij_min_pipe+1), int(ij_max_pipe+1),
                                int(k_min_pipe+1), int(k_max_pipe+1),
                                int(k_min_oph+1), int(k_max_oph+1),
                                oh_perm,
                                LGR_NAME,
                                O)

    # 2. open hole sections
    for idx in drilling_df.index:

        # bbox for open hole section
        oph_columns = ['diameter_m', 'ij_min', 'ij_max', 'k_min', 'k_max', 'oh_perm']
        ID_oph, ij_min_oph, ij_max_oph, k_min_oph, k_max_oph, oh_perm = drilling_df.loc[idx, oph_columns]

        # skipping drilling sections that contain casings
        if k_max_oph < k_max_oph_saved:
            continue

        # 1.1 CARFIN output
        CARFIN_oph(ID_oph,
                    int(ij_min_oph+1), int(ij_max_oph+1),
                    int(ij_min_oph+1), int(ij_max_oph+1),
                    int(k_min_oph+1), int(k_max_oph+1),
                    oh_perm,
                    LGR_NAME,
                    O)

    # 3. cement-bond sections
    k_max_saved = -1
    for idx in casings_df.index:

        # bbox for cement bond
        cb_columns = ['diameter_m', 'ij_min', 'ij_max','toc_k_min', 'toc_k_max', 'cb_perm']
        ID_pipe, ij_min_pipe, ij_max_pipe, toc_k_min_cb, toc_k_max_cb, cb_perm = casings_df.loc[idx, cb_columns]

        # locate the interval of cement-bond within drilling
        for index, row in drilling_df.iterrows():
            kmin, kmax = row['k_min'], row['k_max']

            if toc_k_min_cb >= kmin and toc_k_min_cb < kmax:
                toc_k_min = toc_k_min_cb

            if toc_k_max_cb >= kmin and toc_k_max_cb < kmax:
                toc_k_max = toc_k_max_cb

        for ic, kmax in enumerate(drilling_k_maxs):
            if toc_k_min_cb >= kmax:
                idx_kmin = 
                break
        # the thickness of cement bond
        x_thickness = int((drilling_df.iloc[idx0, 'ij_max'] - casings_df.loc[idx0, 'ij_max']))

        # 2.2 CARFN output
        CARFIN_cement_bond(ID_pipe, 
                            int(ij_min_pipe+1), int(ij_max_pipe+1),
                            int(ij_min_pipe+1), int(ij_max_pipe+1),
                            int(toc_k_min_cb+1), int(toc_k_max_cb+1),
                            x_bd=x_thickness,
                            y_bd=x_thickness,
                            perm=cb_perm,
                            LGR_NAME=LGR_NAME, 
                            O=O)


def df_to_gap_barrier(barriers_mod_df: pd.DataFrame,
                      LGR_NAME: str,
                      O: TextIO) -> None:
    """ convert barrier dataframe to gap format

        Args:
            barriers_mod_df (pd.DataFrame): dataframe for barriers
            LGR_NAME (str): LGR name
            O (TextIO): opened file handle
    """

    for idx, row in barriers_mod_df.iterrows():

        # bbox for barrier
        columns = ['diameter_m', 'ij_min', 'ij_max', 'k_min', 'k_max', 'barrier_perm']
        ID, ij_min_bar, ij_max_bar, k_min_bar, k_max_bar, bar_perm = barriers_mod_df.loc[idx, columns]

        # qc it
        # print('barrier===>', idx, ID, strt_depth, end_depth, barrier_perms)
    
        # 1.1 CARFIN output
        CARFIN_barrier(ID, 
                        int(ij_min_bar+1), int(ij_max_bar+1),
                        int(ij_min_bar+1), int(ij_max_bar+1),
                        int(k_min_bar+1), int(k_max_bar+1),
                        bar_perm,
                        LGR_NAME, 
                        O)



