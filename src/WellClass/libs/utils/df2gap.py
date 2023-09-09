
import pandas as pd

from src.GaP.libs.models import (
    DepthModel,
    PipeCementModel,
    ElemModel,
)

def to_gap_casing_list(drilling_df: pd.DataFrame, 
                        casings_df: pd.DataFrame, 
                        oh_perm: float|int, 
                        cb_perm: float|int) -> list[PipeCementModel]:
    """ convert casing dataframe to gap format

        Args:
            drilling_df (pd.DataFrame): dataframe for drilling
            casings_df (pd.DataFrame): dataframe for casings
            oh_perm (float): permeability for open hole
            cb_perm (float): permeability for cement bond

        Returns:
            list[PipeCementModel]: list of PipeCementModels
    """

    casings_list = []
    for idx, row in casings_df.iterrows():

        ID = row['diameter_m'] 
        strt_depth, end_depth = row['top_msl'], row['bottom_msl']
        strt_depth_cement, end_depth_cement = row['toc_msl'], row['boc_msl']
        strt_depth_oph, end_depth_oph = drilling_df.loc[idx, ['top_msl', 'bottom_msl']]

        # qc it
        # print(idx, ID, strt_depth, end_depth, strt_depth_cement, end_depth_cement, strt_depth_oph, end_depth_oph)
        
        casing_geom = PipeCementModel(ID=ID, 
                                      pipe=DepthModel(strt_depth=strt_depth, end_depth=end_depth, perm=oh_perm),
                                      oph=DepthModel(strt_depth=strt_depth_oph, end_depth=end_depth_oph), 
                                      cement=DepthModel(strt_depth=strt_depth_cement, end_depth=end_depth_cement, perm=cb_perm)
                                    )

        casings_list.append(casing_geom)

    return casings_list

def to_gap_barrier_list(barriers_mod_df: pd.DataFrame, 
                        barrier_perms: list[float]) -> list[ElemModel]:
    """ convert barrier dataframe to gap format
    """

    ib = 0
    barrier_list = []
    for idx, row in barriers_mod_df.iterrows():

        ID, strt_depth, end_depth = row['diameter_m'], row['top_msl'], row['bottom_msl']

        # qc it
        # print(idx, ID, strt_depth, end_depth)
        
        barrier = ElemModel(ID=ID, 
                            pipe=DepthModel(strt_depth=strt_depth, end_depth=end_depth, perm=barrier_perms[ib])
                        )
        ib = ib + 1
        #
        barrier_list.append(barrier)

    return barrier_list

