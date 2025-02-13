"""  
Recipe for generating the data necessary to make the legacy well pressure plots.

Well class examples. In the following, assume they run from project's folder.

--config-file or -c 
  yaml file describing the minimum input data necessary to run the code


        apiVersion: well/v0.1
        kind: Well
        metadata:
        namespace: screen
        name: wildcat
        author: equinor

        spec:

        well_header:

            well_name: 'well_name'
            well_rkb: 99
            sf_depth_msl: 999
            well_td_rkb: 9999
            sf_temp: 9
            geo_tgrad: 99

        
        reservoir_pressure:

            depth_msl: 9999
            RP2: '+ 99'

        co2_datum: 9999

--mixture-index or -mi
  optional value. index to identify the predefined pvt tables, defaulted to 0 = pure CO2

  
--max-pressure-depths or -maxP
  option for depths to eavluate at which CO2 gradient reaches Shmin. evaluation of depth values at which max P is computed.
  it can be a single value, or a list of values.

--out-path or -p
  optional. If the figure wants to be stored, output folder name

--out-name or -o
  optional. If the figure wants to be stored, output file name

--display
  display the figure.


# 1. smeaheia_v1

$ python -m experiments.well_pressure_tables --config-file ./test_data/examples/smeaheia_v1/smeaheia.yaml

# 2. smeaheia_v2

$ python -m experiments.well_pressure_tables --config-file ./test_data/examples/wildcat/wildcat.yaml

"""
import os
import json

import argparse
from argparse import Namespace

from pathlib import Path
import sys

from importlib.resources import files



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
    plot_pressure
    )




def main(args: Namespace):

    ############# 0. User options ######################

    # input configuration file name,
    # for example, './test_data/examples/smeaheia_v1.yaml'
    well_name = Path(args.config_file)

    # input mixture index.
    mixture_index = args.mixture_index
    print(f'{mixture_index=}')

    mixtures_list = ['pure_co2', 'mixture1', 'mixture2']

    mixture = mixtures_list[mixture_index]

    pvt_files = files('src.WellClass.libs.pvt.pvt_constants')
    pvt_path = pvt_files / mixture


    #depth values to evaluyate pressure at Shmin
    maxP_values = args.max_pressure_depths

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
                    co2_datum    = well_csv['co2_datum'],
            )
    
    # 3.2 build pressure class
    my_pressure = Pressure( header      = well_csv['well_header'],
                            reservoir_P = well_csv['reservoir_pressure'],
                            co2_datum   = well_csv['co2_datum'],
                            max_pressure_pos = maxP_values,
                            pvt_path    = pvt_path,)


    # well sketch
    output = None
    if args.out_name:
        if not os.path.exists(args.out_path):
            os.makedirs(args.out_path, exist_ok=True)
        # output file
        output = os.path.join(args.out_path, args.out_name)

    #Plot sketch, pressures
    # plot_sketch_pressure(my_well, my_pressure, save_file=output)
    plot_pressure(my_pressure)

    #Export tables to json files
    my_pressure.pressure_CO2.to_json('pressure_co2.json', orient='records')

    with open('pressure_scenarios.json', 'w') as f:
        # Use json.dump() to write dictionary to a file
        json.dump(my_pressure.pressure_scenarios, f)

    # for qc
    if args.display:
        plt.show()

if __name__ == '__main__':

    # Create the parser
    parser = argparse.ArgumentParser()
    
    #yaml configuration file
    parser.add_argument('-c', "--config-file", type=str, required=True,
                        help="well configuration file, can be .yaml or .csv format")
    

   
    #mixture index, defaulted to 0 = pure CO2
    parser.add_argument('-mi', "--mixture-index", type=int, default=0,
                        help="option of pre-defined chemical mixtures")
    
    #evaluation of depth values at which max P is computed.
    parser.add_argument('-maxP', "--max-pressure-depths", nargs='+', type=float,
                        help="option for depths to eavluate at which CO2 gradient reaches Shmin")


    # save figure to the disk
    parser.add_argument('-p', '--out-path', type=str,
                        default='.',
                        help="output folder for the figure")
    
    parser.add_argument('-o', '--out-name', type=str,
                        default="",
                        help="output file name for the figure")

    # display the figure
    parser.add_argument('--display', action='store_true',
                        help='display of the well sketch')

    # Parse the argument
    args = parser.parse_args()

    main(args)