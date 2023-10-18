#!/usr/bin/env python
#GaP_17062022.py + SATNUM =2 is added for linear relperms inside the pipe (pipe_builder function) 
#GaP_17062022.py + well sketch visualization is added for QC 
#GaP_02052022.py  + the ACTNUM inside CARFIN is changed to TRANX/Y/Z for sake of runtime and visualization 
#GaP_05042022_FPT.py + the grids around legacy well inside LGR are de-activated 
#GaP.py + the reservoir is isolated from OVB 
# the only communication path within OVB and reservoir is pipe (ACTNUM inside CARFIN) 

""" GaP_HIRES_15112022.py - Grid as Pipe is used to generate LGR grid for later CO2 simulation in pflotran or Eclipse

SMEAHEIA case

    script location: /scratch/SCS/eim/SMEAHEIA/GaP_HIRES_15112022.py
    Input: This case takes the files associated to the Pflotran simulation grid called GEN_NOLGR_PH2. It reads its grid (.grdecl file).
    Output: Produces a schematic plot (just informative). And the following file:
    /scratch/SCS/eim/SMEAHEIA/LEG_HIRES.grdecl

    That file is then declared as an include file in the main grdecl file: /scratch/SCS/eim/SMEAHEIA/res_ovb_kdome_nolgr_twoplugs_ph2.grdecl 
    Which is used by the following simulation case: /scratch/SCS/eim/SMEAHEIA/GEN_NOLGR_PH2.in
    You can explore the results of these simulation case using ResInsight.

"""

# To run it, type the following at GaP directory:
#
# 1. To test two barriers,
#
# $ python -m experiments.GaP_HIRES --config-file data/smeaheia/config.yaml --display
#
# 2. To test one barrier,
#
# $ python -m experiments.GaP_HIRES --config-file data/smeaheia_onepipe/config.yaml --display
#

import os
import json
import yaml
import numpy as np

from matplotlib import pyplot as plt
#import rips 

from ecl.grid import EclGrid
from ecl.eclfile import EclInitFile

# 
from libs.geometry import (
    make_stat_3d_main,
    generate_LGR_grids,
)

from libs.carfin import build_grdecl

from libs.plotting import (
    pipe_plotter, 
    barrier_plotter,
    ob_reservoir_plotter,
)

from libs.models import (
    GaPModel,
)

def main(args):

    ###### 0. user-input parameters #########

    # load configuration file
    with open(args.config_file, "r") as f:
        config = yaml.safe_load(f)

    # yaml_config = yaml.dump(config, sort_keys=False)
    # print(yaml_config)

    # dict => pydantic model
    gap = GaPModel(**config)

    # # pretty print
    # print(json.dumps(gap.model_dump(mode='json'), indent=2))

    # alias, only for convenience
    spec = gap.spec

    # 0.1 before sorting
    # casings
    casings = spec.casings
    # barriers
    barriers = spec.barriers

    # # for qc
    # print(casings)

    # 0.2 after sorting
    # casings
    casings.sort(reverse=True, key=lambda elem: elem.ID)
    # barriers
    barriers.sort(reverse=True, key=lambda elem: elem.ID)

    # # for qc
    # print(casings)

    # input simulation case
    input_folder = spec.sim_case.folder
    simcase = os.path.join(input_folder, spec.sim_case.filename)

    # output LGR name
    output_folder = spec.lgr_out.folder
    LGR_NAME = spec.lgr_out.filename

    # minimum dz value for overburden
    MIN_DZ_OB = spec.defaults.mindz_ob


    ###############  1. Initialization ###############

    ####### 1.1 Loading the model #######

    # Eclipse grid
    grid = EclGrid(simcase + ".EGRID") 
    init = EclInitFile(grid, simcase + ".INIT")
    #init = EclGrid(simcase + ".INIT")

    ####### 1.2 Main grid info #######

    # coarse grid
    nx, ny, nz, _ = grid.get_dims()

    # legacy well grid (coarse grid)
    main_grd_i = nx//2
    main_grd_j = ny//2

    main_grd_min_k = 1
    main_grd_max_k = nz

    # dx and dy (coarse grid)
    main_grd_dx = make_stat_3d_main(grid, init, 'DX', 0)[main_grd_i-1, main_grd_j-1, main_grd_min_k-1]       # 656.168 ft 
    main_grd_dy = make_stat_3d_main(grid, init, 'DY', 0)[main_grd_i-1, main_grd_j-1, main_grd_min_k-1]       # 656.168 ft   # noqa: F841

    # list of DZ (coarse grid)
    main_DZ = make_stat_3d_main(grid, init, 'DZ', 0)[main_grd_i-1, main_grd_j-1, (main_grd_min_k-1):main_grd_max_k]

    # # reference depth where LGR (legacy well) starts
    # ref_depth = make_stat_3d_main(grid, init, 'DEPTH', 0)[main_grd_i-1, main_grd_j-1, main_grd_min_k-1] 
    ref_depth = make_stat_3d_main(grid, init, 'DEPTH', 0)[main_grd_i-1, main_grd_j-1, (main_grd_min_k-1)] - \
        (0.5*make_stat_3d_main(grid, init, 'DZ', 0)[main_grd_i-1, main_grd_j-1, (main_grd_min_k-1)])

    # assume any layers that has DZ > MIN_DZ_OB (default: 10) will be overburden layers
    # calculate the number of layerd in overburden
    no_of_layers_in_OB = np.count_nonzero(main_DZ[main_DZ > MIN_DZ_OB])

    ############### 3. generate LGR grids ###############

    print ('Generate LGR grid...') 

    # extract casing IDs
    casing_IDs = list(map(lambda casing: casing.ID, casings))

    # generate LGR grids
    LGR_sizes_xy, LGR_depths, LGR_numb_z, min_grd_size = \
            generate_LGR_grids(casing_IDs,
                                main_grd_dx, main_DZ,
                                ref_depth,
                                main_grd_min_k, main_grd_max_k,
                                no_of_layers_in_OB)

    ############### 4. generate grdecl ###############

    print ('Generate LGR and output...') 

    build_grdecl(output_folder, LGR_NAME,
                    casings,
                    barriers,
                    LGR_sizes_xy, 
                    LGR_depths,
                    LGR_numb_z, 
                    min_grd_size,
                    grid.getNX(), grid.getNY(),
                    main_grd_i, main_grd_j,
                    main_grd_min_k, main_grd_max_k,
                    no_of_layers_in_OB)

    ############### 5. visualize Well Sketch Design ###############
    
    if args.display:

        print ('Visualizes the well and its locationin OVB and reservoir...') 

        # loop for casings
        for casing in casings:
            pipe_plotter (casing,
                            LGR_depths,
                            plt)        

        # 4. barriers
        for barrier in barriers:
            barrier_plotter(barrier,
                            plt)  

        # ref depths for ob and reservoir
        ref_depth1 = make_stat_3d_main(grid, init, 'DEPTH', 0)[main_grd_i-1,main_grd_j-1, 0]
        ref_depth2 = make_stat_3d_main(grid, init, 'DEPTH', 0)[main_grd_i-1,main_grd_j-1, no_of_layers_in_OB]
        # TODO(hzh): 59?
        ref_depth3 = make_stat_3d_main(grid, init, 'DEPTH', 0)[main_grd_i-1,main_grd_j-1, 59]

        # 5. Overburden / Reservoir Design
        ob_reservoir_plotter(ref_depth1, ref_depth2, ref_depth3, plt)

        plt.xlim(-4,4)
        plt.show()  

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--config-file', type=str,
                        default='data/smeaheia/config.yaml',
                        help="yaml configuration file name")
    
    parser.add_argument('--display', action='store_true',
                        help="well sketch when it is on")
    
    args = parser.parse_args()


    main(args)
