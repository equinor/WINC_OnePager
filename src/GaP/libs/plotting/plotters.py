""" well sketch
"""

import numpy as np

from ..models import (
    PipeCementModel,
    ElemModel
)

from ..geometry.bbox_utils import BBoxZ

def pipe_plotter (cg: PipeCementModel,
                  LGR_depths,
                  plt):
    """ casing, cement bonds, and open hole sections
    """

    # for convenience
    ID = cg.ID

    # 1. for Casing lines 
    plt.vlines(x=-ID/2, ymin=-cg.pipe.end_depth, ymax=-cg.pipe.strt_depth) 
    plt.vlines(x=(-ID/2) + ID, ymin=-cg.pipe.end_depth, ymax=-cg.pipe.strt_depth)

    plt.scatter(-ID/2, -cg.pipe.end_depth, marker=8, color='black')
    plt.scatter( ID/2, -cg.pipe.end_depth, marker=9, color='black')

    # for annotations

    # min/max z, pipe
    k_min_pipe, k_max_pipe = BBoxZ(cg.pipe.strt_depth, cg.pipe.end_depth, LGR_depths)

    # min/max z, hole
    k_min_hole, k_max_hole = BBoxZ(cg.oph.strt_depth, cg.oph.end_depth, LGR_depths)
    
    # annotations
    plt.annotate( str(k_max_hole), (-ID/2-0.5, -cg.oph.end_depth))
    plt.annotate( str(k_min_hole), (-ID/2-0.5, -cg.oph.strt_depth))

    plt.annotate( str(k_max_pipe), (ID/2+0.5, -cg.pipe.end_depth))
    plt.annotate( str(k_min_pipe), (ID/2+0.5, -cg.pipe.strt_depth))
    plt.annotate( 'Openhole K index', (-3.5, 5))
    plt.annotate( 'Casing K index', (+2,0))

    # 2. for cement bonds

    # 2.1 cement bond thickness is hardcoded to 0.05 for visualization 
    plt.fill_between ([-ID/2, -ID/2 - 0.08], 
                      -cg.cement.end_depth,
                      -cg.cement.strt_depth, 
                      color='grey', 
                      alpha = 0.5)
    
    # 2.2 cement bond thickness is hardcoded to 0.05 for visualization 
    plt.fill_between ([ID/2, ID/2 + 0.08 ], 
                      -cg.cement.end_depth,
                      -cg.cement.strt_depth, 
                      color='grey', 
                      alpha = 0.5)
    
    # 3. for open hole section 
    plt.fill_between([-ID/2, (-ID/2)+ID], 
                     -cg.oph.end_depth,
                     -cg.oph.strt_depth, 
                     color = 'green', 
                     alpha = 0.2)

def barrier_plotter (cg: ElemModel,
                     plt):
    """ barrier
    """

    # for convenience
    ID = cg.ID

    # 4. for barrier
    plt.fill_between(  [-ID/2, (-ID/2) + ID],
                        -cg.pipe.end_depth , 
                        -cg.pipe.strt_depth, 
                        color = 'red', 
                        alpha = 0.5)

def ob_reservoir_plotter(ref_depth1, ref_depth2, ref_depth3, plt):
    """ overburden and reservoir
    """
    # 5.1 OB
    plt.fill_between([-4,4], -ref_depth1, -ref_depth2, color='blue', alpha=0.1)

    # 5.2 reservoir
    plt.fill_between([-4,4], -ref_depth2, -ref_depth3, color='brown', alpha=0.5)

