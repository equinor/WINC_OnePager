
import os
import sys
import json

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt



# where WellClass and Ga[ codes are located
sys.path.append('../')


# GaP
from src.GaP.libs.models import (
    PipeCementModel,
    ElemModel,
    DepthModel,
)

from src.GaP.libs.carfin import build_grdecl

# WellClass
from src.WellClass.libs.well_class import Well

from src.WellClass.libs.utils import (
    csv_parser,
    yaml_parser,
)
from src.WellClass.libs.utils.LGR_grid import (
    compute_ngrd,
    generate_LGR_xy,
    generate_LGR_z,
)

from src.WellClass.libs.utils.df2gap import (
    gap_casing_list,
    gap_barrier_list
)


# plots
from src.WellClass.libs.plotting import plot_well_perm

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
               'sim_path': r'../test_data/examples/smeaheia_v1',
               'simcase': r'GEN_NOLGR_PH2'}
smeaheia_v2 = {'well_input': r'GaP_input_Smeaheia_v3.csv', 
               'well_input_yaml': r'smeaheia.yaml', 
            #    'sim_path': r'/scratch/SCS/bkh/wbook/realization-0/iter-0/pflotran/model', 
               'sim_path': r'../test_data/examples/smeaheia_v2', 
               'simcase': r'TEMP-0'}
cosmo = {
         'well_input': r'GaP_input_Cosmo_v3.csv', 
         'well_input_yaml': r'cosmo.yaml', 
        #  'sim_path': r'/scratch/SCS/bkh/well_class_test1/realization-0/iter-0/pflotran/model', 
         'sim_path': r'../test_data/examples/cosmo', 
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



############################# grid_init ################




##############
# # LGR grid information in x, y, z directions
# 
# We are going to compute the grid sizes in lateral (x and y) and vertical directions

# ### 1. Compute minimum grid size

# 0. minimum grid size

# minimum grid size depends on minimum annulus thickness
min_grd_size = casings_df['thick_m'].min()

if min_grd_size < 0.05:
    min_grd_size = 0.05

print(f'Minimimum grid size is {min_grd_size*100:.2f} cm')

# TODO(hzh): manually set it
if Ali_way:
    min_grd_size = 0.05

# 1. compute number of LGR grids for drilling, casing and borehole, respectively

# for drilling
drilling_df['n_grd_id']  = drilling_df['diameter_m'].map(lambda x: compute_ngrd(x, min_grd_size))

# the following two are not used
casings_df[ 'n_grd_id']  = casings_df['diameter_m'].map(lambda x: compute_ngrd(x, min_grd_size))
borehole_df['n_grd_id'] = borehole_df['id_m'].map(lambda x: compute_ngrd(x, min_grd_size))

# ### 2. Compute LGR grid sizes in x-y directions

# 2.1 Number of cells representing horizontal LGR
no_latral_fine_grd = drilling_df['n_grd_id'].max()

# 2.2 generate the LGR grid sizes in x-y
LGR_sizes_x, LGR_sizes_y, _ = generate_LGR_xy(no_latral_fine_grd, 
                                              min_grd_size, 
                                              main_grd_dx, main_grd_dy,
                                              Ali_way=Ali_way)

# ### 3. Compute LGR grid sizes in z direction

# LGR_sizes_z, LGR_numb_z, LGR_depths, _ = generate_LGR_z(DZ_rsrv, DZ_ovb_coarse)
# TODO(hzh): to make LGR starts at ref_depth
LGR_sizes_z, LGR_numb_z, LGR_depths, _ = generate_LGR_z(DZ_rsrv, DZ_ovb_coarse, ref_depth)


#####################################
# # Set up dataframe for LGR mesh

########################################33333
# # Bounding box for well elements

# # Write LGR file

output_dir = '.'

# LRG name 
LGR_NAME = 'LEG_HIRES'

# prepare info about Casing, Cement Bond and Open hole  for GaP
casing_list = gap_casing_list(drilling_df, 
                              casings_df, 
                              oh_perm, cb_perm)

# prepare info about Barrier for GaP 
barrier_list = gap_barrier_list(barriers_mod_df, 
                                barrier_perms)

# generate .grdecl file
# TODO(hzh): add 1s to indices here
build_grdecl(output_dir, LGR_NAME,
                casing_list,
                barrier_list,
                LGR_sizes_x, 
                LGR_depths, 
                LGR_numb_z, 
                min_grd_size,
                grid.getNX(), grid.getNY(),
                main_grd_i + 1, main_grd_j + 1,
                main_grd_min_k + 1, main_grd_max_k + 1,
                no_of_layers_in_OB)



