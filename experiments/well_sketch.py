

"""  Well sketch examples,

# 1. smeaheia_v1

$ python -m experiments.well_sketch --sim-path ./test_data/examples/smeaheia_v1 --well smeaheia.yaml --plot 

# 2. smeaheia_v2

$ python -m experiments.well_sketch --sim-path ./test_data/examples/smeaheia_v2 --well smeaheia.yaml --plot

# 3. wildcat

$ python -m experiments.well_sketch --sim-path ./test_data/examples/wildcat --well wildcat.yaml --plot

"""

import json

import argparse
from argparse import Namespace

import pathlib

from matplotlib import pyplot as plt

from src.WellClass.libs.utils import (
    csv_parser,
    yaml_parser,
)

# WellClass
from src.WellClass.libs.well_class import Well

# plot
from src.WellClass.libs.plotting import (
    plot_sketch
)

def main(args: Namespace):

    ############# 0. User options ######################

    # where the location for the input parameters and eclipse .EGRID and .INIT files
    # configuration path, for example './test_data/examples/smeaheia_v1'
    sim_path = pathlib.Path(args.sim_path)

    # input configuration file name, for example, 'smeaheia.yaml'
    well_config = pathlib.Path(args.well)

    # extract suffix
    suffix = well_config.suffix
    # .yaml or .csv?
    use_yaml = suffix in ['.yaml', '.yml']

    ############ 2. Load well configuration file ###############

    # where well configuration file is located
    well_name = sim_path/well_config
    
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
    
    # well sketch
    plot_sketch(my_well, save_file=args.file_name)

    # for qc
    if args.plot:

        plt.show()

if __name__ == '__main__':

    # Create the parser
    parser = argparse.ArgumentParser()

    parser.add_argument('-p', "--sim-path", type=str, required=True,
                        help='The file path to the configuration folder')
    
    parser.add_argument('-w', "--well", type=str, required=True,
                        help="input well configuration file name, can be .yaml or .csv format")
    
    parser.add_argument('--plot', action='store_true', help='plot well sketch and well grids')

    # save figure to the disk
    parser.add_argument('-o', '--file_name', type=str, help="output file name for the figure")

    # Parse the argument
    args = parser.parse_args()

    main(args)

