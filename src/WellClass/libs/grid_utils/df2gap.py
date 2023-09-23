
from typing import List, Union

import pandas as pd

from src.GaP.libs.models import (
    DepthModel,
    PipeCementModel,
    ElemModel,
)

def to_gap_casing_list(drilling_df: pd.DataFrame, 
                        casings_df: pd.DataFrame) -> List[Union[PipeCementModel, ElemModel]]:
    """ convert casing dataframe to gap format

        Args:
            drilling_df (pd.DataFrame): dataframe for drilling
            casings_df (pd.DataFrame): dataframe for casings

        Returns:
            list[PipeCementModel]: list of PipeCementModels
    """

    casings_list = []
    for idx, row in drilling_df.iterrows():

        ID_oph, strt_depth_oph, end_depth_oph, oh_perm = drilling_df.loc[idx, ['diameter_m', 'top_msl', 'bottom_msl', 'oh_perm']]

        try:
            # casings
            ID, strt_depth, end_depth = casings_df.loc[idx, ['diameter_m', 'top_msl', 'bottom_msl']]
            # cement
            strt_depth_cement, end_depth_cement, cb_perm = casings_df.loc[idx, ['toc_msl', 'boc_msl', 'cb_perm']]
        except:
            print('DEBUG ==============> handling open hole section')
            casing_geom = ElemModel(ID=ID_oph,
                                    pipe=DepthModel(strt_depth=strt_depth_oph, end_depth=end_depth_oph, perm=oh_perm))
            casings_list.append(casing_geom)

            continue

        # qc it
        # print('casing===>', idx, ID, strt_depth, end_depth, strt_depth_cement, end_depth_cement, strt_depth_oph, end_depth_oph)
        # print('casing===>', idx, ID, ', cb_perm=', cb_perm)
        
        casing_geom = PipeCementModel(ID=ID, 
                                      pipe=DepthModel(strt_depth=strt_depth, end_depth=end_depth, perm=oh_perm),
                                      oph=DepthModel(strt_depth=strt_depth_oph, end_depth=end_depth_oph), 
                                      cement=DepthModel(strt_depth=strt_depth_cement, end_depth=end_depth_cement, perm=cb_perm)
                                    )

        casings_list.append(casing_geom)

    return casings_list

def to_gap_barrier_list(barriers_mod_df: pd.DataFrame) -> List[ElemModel]:
    """ convert barrier dataframe to gap format
    """

    ib = 0
    barrier_list = []
    for idx, row in barriers_mod_df.iterrows():

        ID, strt_depth, end_depth, barrier_perms = row['diameter_m'], row['top_msl'], row['bottom_msl'], row['barrier_perm']

        # qc it
        # print('barrier===>', idx, ID, strt_depth, end_depth, barrier_perms)
        
        barrier = ElemModel(ID=ID, 
                            pipe=DepthModel(strt_depth=strt_depth, end_depth=end_depth, perm=barrier_perms)
                        )
        ib = ib + 1
        #
        barrier_list.append(barrier)

    return barrier_list

