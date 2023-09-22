
""" This module will generate .grdecl file for lgr grid. It uses pflotran's dry runs to generate coarse grid information and lgr grid information.

$ python -m experiments.gap_pflotran \
    --sim-path ./test_data/examples/cosmo-pflotran \
    --well cosmo.yaml \
    --sim-case1 TEMP-0_NOSIM \
    --sim-case2 TEMP-0 \
    --plot

"""

import os
import json

import argparse
import pathlib
import subprocess

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
    GridLGR,
    LGRBuilder,
)

# plots
from src.WellClass.libs.plotting import (
    plot_grid,
)

def main(args):

    ############# 0. User options ######################

    # TODO(hzh): use Ali's grid logic
    Ali_way = args.ali_way

    # where the location for the input parameters and eclipse .EGRID and .INIT files
    # configuration path, for example './test_data/examples/cosmo-pflotran'
    sim_path = pathlib.Path(args.sim_path)

    # configuration filename
    well_input = pathlib.Path(args.well)

    # dry run: filename prefix on coarse grid, e.g., TEMP-0_NOSIM
    sim_case_NOSIM = args.sim_case1

    # lgr run: filename prefix on lgr grid, e.g., TEMP-0
    sim_case_LGR = args.sim_case2

    ############# 1. computed parameters ######################

    # extract suffix from the configuration file name
    file_extension = well_input.suffix

    # .yaml or .csv?
    use_yaml = False
    if file_extension in ['.yaml', '.yml']:
        use_yaml = True

    # file prefix for dry run
    # where eclipse .EGRID and .INIT files will be located
    simcase1 = sim_path/'model'/sim_case_NOSIM
    simcase1

    # LGR
    simcase2 = sim_path/'model'/sim_case_LGR
    simcase2

    ############ 1.5 Run coarse simulation ######################
    # file name (coarse grid) for pflotran run
    run_config_coarse = simcase1.with_suffix('.in')

    # the command
    run_command = f'runpflotran1.8 -i -nm 6 {run_config_coarse}'
    command_array = run_command.split()

    # launch the command
    results = subprocess.run(command_array, capture_output=True, encoding="utf-8")
    print(results.stdout)

    ############ 2. Load well configuration file ###############

    # where well configuration file is located
    well_name = sim_path/well_input
    
    if use_yaml:
        
        # # pydantic model
        well_model = yaml_parser(well_name)
        well_csv = json.loads(well_model.spec.model_dump_json())
    else:

        # load the well information
        well_csv = csv_parser(well_name)

    ########### 3. build Well class ######################

    # 3.1 build well class
    my_well = Well( header       = well_csv['well_header'], 
                    drilling     = well_csv['drilling'],
                    casings      = well_csv['casing_cement'],
                    geology      = well_csv['geology'],
                    barriers     = well_csv['barriers'], 
                    barrier_perm = well_csv['barrier_permeability'],
                    co2_datum    = well_csv['co2_datum'],
            )

    # 3.2 to well dataframe
    well_df = WellDataFrame(my_well)

    # for convenience

    # 3.3 extract dataframes
    annulus_df = well_df.annulus_df
    drilling_df = well_df.drilling_df
    casings_df = well_df.casings_df
    borehole_df = well_df.borehole_df

    barriers_mod_df = well_df.barriers_mod_df

    ############### 4. various grids #####################

    ##### 4.1 grid_coarse 

    # Loading the model from .EGRID and .INIT
    grid_coarse = GridCoarse(str(simcase1))

    ##### 4.2 LGR grid 

    # LGR grid information in x, y, z directions
    lgr = LGRBuilder(grid_coarse, 
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

    # set up LGR grid
    grid_refine.build_LGR(drilling_df, casings_df, barriers_mod_df)


    ########### 6. output grdecl file ###################

    LGR_NAME = 'TEMP_LGR'
    output_dir = sim_path/'include'

    # Write LGR file
    lgr.build_grdecl(output_dir, LGR_NAME,
                     drilling_df,
                     casings_df,
                     barriers_mod_df)
    
    ########### 7. Run LGR simulation ###################
    # file name (LGR grid) for pflotran run
    run_config_lgr = simcase2.with_suffix('.in')

    # the command
    run_command = f'runpflotran1.8 -i -nm 6 {run_config_lgr}'
    command_array = run_command.split()

    # launch the command
    results = subprocess.run(command_array, capture_output=True, encoding="utf-8")
    print(results.stdout)

    # for qc
    if args.plot:

        # load LGR grid from simulation file
        grid_lgr = GridLGR(str(simcase2))

        # coarse grid
        plot_grid(my_well, grid_coarse)

        # LGR grid from dataframe
        plot_grid(my_well, grid_refine)

        # LGR grid from pflotran output
        plot_grid(my_well, grid_lgr)

if __name__ == '__main__':

    # Create the parser
    parser = argparse.ArgumentParser()

    parser.add_argument("--ali-way", action="store_true",
                        help="Use Ali's logic to generate LGR grids")

    parser.add_argument('-p', "--sim-path", type=str, required=True,
                        help='The file path to the configuration folder')
    
    parser.add_argument('-w', "--well", type=str, required=True,
                        help="input well configuration file name, can be .yaml or .csv format")

    parser.add_argument('-s1', "--sim-case1", type=str, required=True,
                        help="file name prefix for dry run")

    parser.add_argument('-s2', "--sim-case2", type=str, required=True,
                        help="file name prefix  for lrg run")
    
    parser.add_argument('--plot', action='store_true', help='plot well sketch and well grids')

    # Parse the argument
    args = parser.parse_args()

    main(args)

