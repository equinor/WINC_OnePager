

from src.GaP.libs.geometry.bbox_pipe_with_oph_sec import bbox_for_pipe_with_oph_sec
from src.GaP.libs.geometry.bbox_cement_bond import bbox_for_cement_bond
from src.GaP.libs.geometry.bbox_barrier import bbox_for_barrier

def test_casing_bbox(casing_fixture, 
                     LGR_fixture, 
                     bbox_fixture,
                     ):
    """ test bbox for casing
    """
    casing_geom = casing_fixture
    LGR_sizes_xy, LGR_depths, min_grd_size = LGR_fixture
    bbox_casing = bbox_fixture[0]

    bbox_casing_test = bbox_for_pipe_with_oph_sec (casing_geom, 
                                                    LGR_sizes_xy, 
                                                    LGR_depths, 
                                                    min_grd_size)
    
    for i in range(len(bbox_casing)):
        assert bbox_casing[i] == bbox_casing_test[i]

def test_cement_bbox(casing_fixture, 
                     LGR_fixture, 
                     bbox_fixture,
                     ):
    """ test bbox for cement bond in the casing
    """
    casing_geom = casing_fixture
    LGR_sizes_xy, LGR_depths, min_grd_size = LGR_fixture
    bbox_cement_bond = bbox_fixture[1]

    bbox_cement_bond_test = bbox_for_cement_bond (casing_geom, 
                                                    LGR_sizes_xy, 
                                                    LGR_depths, 
                                                    min_grd_size)
    
    for i in range(len(bbox_cement_bond)):
        assert bbox_cement_bond[i] == bbox_cement_bond_test[i]

def test_barrier_bbox(barrier_fixture, 
                      LGR_fixture, 
                      bbox_fixture,
                      ):
    """ test bbox for barrier
    """
    barrier_geom = barrier_fixture
    LGR_sizes_xy, LGR_depths, min_grd_size = LGR_fixture
    bbox_barrier = bbox_fixture[2]

    bbox_barrier_test = bbox_for_barrier (barrier_geom, 
                                            LGR_sizes_xy, 
                                            LGR_depths, 
                                            min_grd_size)
    for i in range(len(bbox_barrier)):
        assert bbox_barrier[i] == bbox_barrier_test[i]