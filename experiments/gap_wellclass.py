
""" To run it, type the following:

$ python -m experiments.gap_wellclass --use-yaml --case-type smeaheia_v1 --plot --ali-way 

for comparing the output from smeaheia data with the output using Ali's grid logic.

Otherwise, use the following:

$ python -m experiments.gap_wellclass --use-yaml --case-type cosmo --plot 

"""
import os
import json

import argparse

from src.WellClass.libs.utils import (
    csv_parser,
    yaml_parser,
)

# WellClass
from src.WellClass.libs.well_class import Well

from src.WellClass.libs.grid_utils import (
    WellDataFrame,
    GridCoarse,
    GridRefine,
    LGR,
)

# plots
from src.WellClass.libs.plotting.plot_grids import (
    plot_coarse,
    plot_refine,
)

# # Examples
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

def main(args):

    ############# 0. User options ######################

    # TODO(hzh): use Ali's grid logic
    Ali_way = args.ali_way

    # use yaml or csv input file
    use_yaml = args.use_yaml

    # pick an example from given three options: 
    #  i.e, smeaheia_v1, smeaheia_v2 and cosmo
    case_type = args.case_type

    # output
    output_dir = args.output_dir
    # LRG name 
    LGR_NAME = args.output_name


    ############# 1. Selected case ####################

    # the selected example for testing
    case = examples[case_type]

    # where the location for the input parameters and eclipse .EGRID and .INIT files
    sim_path = case['sim_path']


    ############ 2. Load well configuration file ###############

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

    ########### 3. build Well class ######################

    # build well class
    my_well = Well( header       = well_csv['well_header'], 
                    drilling     = well_csv['drilling'],
                    casings      = well_csv['casing_cement'],
                    geology      = well_csv['geology'],
                    barriers     = well_csv['barriers'], 
                    barrier_perm = well_csv['barrier_permeability'],
                    co2_datum    = well_csv['co2_datum'],
            )

    # to well dataframe
    well_df = WellDataFrame(my_well)

    # for convenience

    # extract dataframes
    annulus_df = well_df.annulus_df
    drilling_df = well_df.drilling_df
    casings_df = well_df.casings_df
    borehole_df = well_df.borehole_df

    barriers_mod_df = well_df.barriers_mod_df

    ############### 4. various grids #####################

    ##### 4.1 grid_coarse 

    # location of .egrid
    simcase = os.path.join(sim_path, case['simcase'])

    # Loading the model
    grid_coarse = GridCoarse(simcase)

    ##### 4.2 LGR grid 

    # LGR grid information in x, y, z directions
    lgr = LGR(grid_coarse, 
              annulus_df, 
              drilling_df, casings_df, borehole_df,
              Ali_way)

    ##### 4.3 grid refine 
    
    # Set up dataframe for LGR mesh
    grid_refine = GridRefine(grid_coarse,
                            lgr.LGR_sizes_x, lgr.LGR_sizes_y, 
                            lgr.LGR_sizes_z, 
                            )

    ############### 5. build LGR #####################

    ##### 5.1 compute bounding box 

    # Bounding box for well elements
    well_df.compute_bbox(grid_refine.mesh_df, grid_refine.nx)

    ##### 5.2 set up material type 

    # set up material type
    grid_refine.set_material_type(drilling_df,
                                  casings_df, 
                                  barriers_mod_df)
    
    ##### 5.3 set up permeability 

    # set up permeability
    grid_refine.set_permeability(drilling_df, casings_df, barriers_mod_df)

    ########### 6. output grdecl file ###################

    # Write LGR file
    lgr.build_grdecl(output_dir, LGR_NAME,
                     drilling_df,
                     casings_df,
                     barriers_mod_df)
    
    # for qc
    if args.plot:
        plot_coarse(my_well, grid_coarse)
        plot_refine(my_well, grid_refine)

if __name__ == '__main__':

    # Create the parser
    parser = argparse.ArgumentParser()

    parser.add_argument("--ali-way", action="store_true",
                        help="Use Ali's logic to generate LGR grids")

    parser.add_argument("--use-yaml", action="store_true",
                        help="Use yaml format as input configuration file")

    parser.add_argument("--case-type", type=str, default='smeaheia_v1',
                        choices=['smeaheia_v1', 'smeaheia_v2', 'cosmo'],
                        help="name of test example")
            
    parser.add_argument('--output-dir', type=str, default='./experiments',
                        help="output directory")
    parser.add_argument('--output-name', type=str, default='LEG_HIRES', 
                        help='output file name')
    
    parser.add_argument('--plot', action='store_true', help='plot well sketch and well grids')

    # Parse the argument
    args = parser.parse_args()

    main(args)

