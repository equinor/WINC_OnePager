
""" To run it, type the following:

$ python -m experiments.gap_wellclass

"""
import os
import sys
import json

import matplotlib.pyplot as plt

# where WellClass and Ga[ codes are located
# sys.path.append('../')


from src.GaP.libs.carfin import build_grdecl

# WellClass
from src.WellClass.libs.well_class import Well

from src.WellClass.libs.utils import (
    csv_parser,
    yaml_parser,
)

from src.WellClass.libs.utils.df2gap import (
    to_gap_casing_list,
    to_gap_barrier_list
)


# plots
from .plot_utils import (
    plot_coarse,
    plot_refine,
)

# 
from .well_df import WellDataFrame

from .grid_coarse import GridCoarse
from .grid_refine import GridRefine
from .LGR    import LGR

def main():

    # ## Some user options

    # TODO(hzh): use Ali's algorithm
    Ali_way = True

    # use yaml or csv input file
    use_yaml = True

    # pick an example from given three options

    case_type = 'cosmo'

    case_type = 'smeaheia_v1'

    # case_type = 'smeaheia_v2'

    # # Examples
    # 
    # The following are the test examples.

    # examples
    smeaheia_v1 = {'well_input': r'GaP_input_Smeaheia_v3.csv', 
                'well_input_yaml': r'smeaheia.yaml', 
                #    'sim_path': r'/scratch/SCS/eim/SMEAHEIA', 
                'sim_path': r'./test_data/examples/smeaheia_v1',
                'simcase': r'GEN_NOLGR_PH2'}
    smeaheia_v2 = {'well_input': r'GaP_input_Smeaheia_v3.csv', 
                'well_input_yaml': r'smeaheia.yaml', 
                #    'sim_path': r'/scratch/SCS/bkh/wbook/realization-0/iter-0/pflotran/model', 
                'sim_path': r'./test_data/examples/smeaheia_v2', 
                'simcase': r'TEMP-0'}
    cosmo = {
            'well_input': r'GaP_input_Cosmo_v3.csv', 
            'well_input_yaml': r'cosmo.yaml', 
            #  'sim_path': r'/scratch/SCS/bkh/well_class_test1/realization-0/iter-0/pflotran/model', 
            'sim_path': r'./test_data/examples/cosmo', 
            'simcase': r'TEMP-0'}

    examples = {
        'smeaheia_v1': smeaheia_v1,
        'smeaheia_v2': smeaheia_v2,
        'cosmo': cosmo
    }

    # # Load well CSV or yaml configuration file
    # 
    # Process CSV with well class.
    # Predefine a dictionary that includes the input CSV well file, the simulation path, and the PFT sim case name

    # the selected example for testing
    case = examples[case_type]

    # root_path = '/scratch/SCS/gpb/SCREEN/GaP_code'

    # where the location for the input parameters and eclipse .EGRID and .INIT files
    sim_path = case['sim_path']

    if use_yaml:
        # where well configuration file is located
        well_name = os.path.join(sim_path, case['well_input_yaml'])
        
        # # pydantic model
        well_model = yaml_parser(well_name)
        well_csv = json.loads(well_model.spec.model_dump_json())
    else:
        # where well configuration file is located
        well_name = os.path.join(sim_path, case['well_input'])

        # load the well information
        well_csv = csv_parser(well_name)

    #Process well by running well class
    my_well = Well( header       = well_csv['well_header'], 
                    drilling     = well_csv['drilling'],
                    casings      = well_csv['casing_cement'],
                    geology      = well_csv['geology'],
                    barriers     = well_csv['barriers'], 
                    barrier_perm = well_csv['barrier_permeability'],
                    co2_datum    = well_csv['co2_datum'],
            )

    # well dataframe
    well_df = WellDataFrame(my_well)

    # dataframes
    annulus_df = well_df.annulus_df
    drilling_df = well_df.drilling_df
    casings_df = well_df.casings_df
    borehole_df = well_df.borehole_df

    barriers_df = well_df.barriers_df
    barriers_mod_df = well_df.barriers_mod_df

    # # Loading the model
    # 
    # - Load the PFT grid, init and restart files.
    # - Grid contains geometry specs
    # - INIT contains static properties (i.e. poro., perm., transmissibilities)
    # - RST contains dynamic properties (i.e. saturations, pressure)
    # 

    # simulation case without legacy well 

    # path = '/scratch/SCS/bkh/wbook/realization-0/iter-0/pflotran/model'

    # location of .egrid
    simcase = os.path.join(sim_path, case['simcase'])

    ############################# grid_coarse ################
    grid_coarse = GridCoarse(simcase)

    ##############
    # # LGR grid information in x, y, z directions
    # 
    # We are going to compute the grid sizes in lateral (x and y) and vertical directions

    lgr = LGR(grid_coarse, 
              annulus_df, 
              drilling_df, casings_df, borehole_df,
              Ali_way)

    #####################################
    # # Set up dataframe for LGR mesh

    grid_refine = GridRefine(grid_coarse,
                            lgr.LGR_sizes_x, lgr.LGR_sizes_y, 
                            lgr.LGR_sizes_z, 
                            )

    #############################################

    # # Bounding box for well elements
    well_df.compute_bbox(grid_refine.mesh_df, grid_refine.nx)

    #################################################
    # set up material type

    grid_refine.set_material_type(drilling_df,
                                  casings_df, 
                                  barriers_df)
    
    ##############################
    # set up permeability
    

    # open hole
    oh_perm = 1e5

    # cemont bond
    cb_perm = 0.5

    # barrier
    barrier_perm = 0.5

    # TODO(hzh): test with original permeabilities

    # open hole
    oh_perm = 10000
    # cemont bond
    cb_perm = 5

    # barrier
    barrier_perm = 0.5

    # TODO(hzh): need to figure out a better to load barrier perms.
    # here manually set them
    barriers_defaults = [0.5, 100, 5, 5, 10]

    # barriers
    barrier_perms = []
    for i in range(len(barriers_mod_df)):
        barrier_perms.append(barriers_defaults[i])

    # set up permeability
    grid_refine.set_permeability(oh_perm, cb_perm, barrier_perm)

    #############################################

    # # Write LGR file

    output_dir = '.'

    # LRG name 
    LGR_NAME = 'LEG_HIRES'

    # prepare info about Casing, Cement Bond and Open hole  for GaP
    casing_list = to_gap_casing_list(drilling_df, 
                                     casings_df, 
                                     oh_perm, cb_perm)

    # prepare info about Barrier for GaP 
    barrier_list = to_gap_barrier_list(barriers_mod_df, 
                                       barrier_perms)

    # generate .grdecl file
    # TODO(hzh): add 1s to indices here
    build_grdecl(output_dir, LGR_NAME,
                 casing_list,
                 barrier_list,
                 lgr.LGR_sizes_x, 
                 lgr.LGR_depths, 
                 lgr.LGR_numb_z, 
                 lgr.min_grd_size,
                 lgr.NX, lgr.NY,
                 lgr.main_grd_i + 1, lgr.main_grd_j + 1,
                 lgr.main_grd_min_k + 1, lgr.main_grd_max_k + 1,
                 lgr.no_of_layers_in_OB)

    plot_coarse(my_well, grid_coarse)
    plot_refine(my_well, grid_refine)

if __name__ == '__main__':

    main()

