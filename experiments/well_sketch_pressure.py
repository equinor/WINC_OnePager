

"""  Well sketch examples. In the following, assume they run from project's folder.

# 1. smeaheia_v1

$ python -m experiments.well_sketch_pressure --config-file ./test_data/examples/smeaheia_v1/smeaheia.yaml -pvt ./test_data/pvt_constants

# 2. smeaheia_v2

$ python -m experiments.well_sketch_pressure --config-file ./test_data/examples/smeaheia_v2/smeaheia.yaml -pvt ./test_data/pvt_constants

# 3. wildcat

$ python -m experiments.well_sketch_pressure --config-file ./test_data/examples/wildcat/wildcat.yaml -pvt ./test_data/pvt_constants

"""
import os
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

# Pressure Class
from src.WellClass.libs.well_pressure import Pressure


# plotting libraries
from src.WellClass.libs.plotting import (
    plot_sketch_pressure
    )

from src.WellClass.tools import init_well_pressure

def main(args: Namespace):

    ############# 0. User options ######################

    # input configuration file name,
    # for example, './test_data/examples/smeaheia_v1.yaml'
    well_name = pathlib.Path(args.config_file)

    pvt_path = pathlib.Path(args.pvt_db_path)

    # extract suffix
    suffix = well_name.suffix
    # .yaml or .csv?
    use_yaml = suffix in ['.yaml', '.yml']

    ############ 2. Load well configuration file ###############
    
    if use_yaml:
        
        # # pydantic model
        well_model = yaml_parser(well_name)
        well_csv = json.loads(well_model.spec.model_dump_json())
    else:

        # load the well information
        well_csv = csv_parser(well_name)

    ########### 3. build Well and pressure classes ######################

    # 3.1 build well class
    my_well = Well( header       = well_csv['well_header'], 
                    drilling     = well_csv['drilling'],
                    casings      = well_csv['casing_cement'],
                    geology      = well_csv['geology'],
                    barriers     = well_csv['barriers'], 
                    barrier_perm = well_csv['barrier_permeability'],
                    co2_datum    = well_csv['co2_datum'],
            )
    
    # 3.2 build pressure class
    my_pressure = Pressure( header      = well_csv['well_header'],
                            reservoir_P = well_csv['reservoir_pressure'],
                            co2_datum   = well_csv['co2_datum'],
                            pvt_path    = pvt_path,)


    # well sketch
    output = None
    if args.out_name:
        if not os.path.exists(args.out_path):
            os.makedirs(args.out_path, exist_ok=True)
        # output file
        output = os.path.join(args.out_path, args.out_name)

    #Plot sketch, pressures
    plot_sketch_pressure(my_well, my_pressure, save_file=output)


    # for qc
    if not args.nodisplay:

        plt.show()

if __name__ == '__main__':

    # Create the parser
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-c', "--config-file", type=str, required=True,
                        help="well configuration file, can be .yaml or .csv format")
    

    parser.add_argument('-pvt', "--pvt-db-path", type=str, required=True,
                        help="pvt fluid database path")
    
    # save figure to the disk
    parser.add_argument('-p', '--out-path', type=str,
                        default='.',
                        help="output folder for the figure")
    
    parser.add_argument('-o', '--out-name', type=str,
                        default="",
                        help="output file name for the figure")

    # display the figure
    parser.add_argument('--nodisplay', action='store_true',
                        help='no display of the well sketch')

    # Parse the argument
    args = parser.parse_args()

    main(args)

