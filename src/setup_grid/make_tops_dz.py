#!/usr/bin/env python
'''
Creates the TOPS and DZ section for the main grid based on input from design-matrix

The results will be zometing like this:

EQUALS 
TOPS  4  4* 1 1 /
/

DZ 
400*40 400*110 3200*150  20000*5 /

Author: bk hegstad, December 2023

'''


import numpy as np
import fmu.config.utilities as utils
import argparse
import pathlib


#If water depth is more than half of the overburden cell thickness (dz) then don't do the "split" of the first top ob layer.
WATER_DEPTH_OB_CRITERIA = 0.5


parser = argparse.ArgumentParser(description="Calculates the grid dimensions for the grdecl-file")

parser.add_argument("outfile", help="path and file to output data")
args = parser.parse_args()


#Read data from design matrix/global variables
CFG = utils.yaml_load("./fmuconfig/output/global_variables.yml")["global"]["GRID_GEO"]
#print(CFG)

tops        = CFG['TOPS']

layer_noc   = CFG['GRD_NX']*CFG['GRD_NY']

water_depth = CFG['WATER_DEPTH']
water_nz    = CFG['WATER_NZ']

ob_nz       = CFG['OB_NZ']

res_nz      = CFG['RES_NZ']
res_dz      = CFG['RES_DZ']
res_depth   = CFG['RES_DEPTH']


#Make the TOPS part
result  =  "EQUALS\n"
result += f"TOPS {int(tops)} 4* 1 1 /\n/\n\n"

#Normalize, only interested in delta-depths from tops.
res_depth   -= tops

#Getting overburden dz. This first ob layer is later split into water and ob
ob_dz = res_depth/ob_nz

#Getting the water part
result += "DZ\n"
result += f"{layer_noc}*{water_depth}"
cum_dz = water_depth

#Getting the first OB layer thickness
dz = ob_dz - water_depth
if dz > WATER_DEPTH_OB_CRITERIA*ob_dz:  #If water layer is too thick - do not do this split 
    result += f" {layer_noc}*{dz:.0f} " #The first layer in the overburden
    ob_nz  -= 1                         #Now the first OB layer is used to adjust to the water depth. Not sure why we would like to do this
    cum_dz += dz

#Getting rest of the overburden
dz = (res_depth - cum_dz)/ob_nz               #Remaining depth interval to fill with cells
result += f" {layer_noc*ob_nz:.0f}*{dz:.0f} "

#Getting the reservoir
result += f" {layer_noc*res_nz:.0f}*{res_dz:.0f} /\n\n\n"



print(result)

#Check if path exists. If not create it
path = pathlib.Path(args.outfile).parent
try:
    path.mkdir(parents=True, exist_ok=False)
except FileExistsError:
    print(f"Folder {str(path)} exists")
else:
    print("Folder {str(path)} is created")

print(f"Writing this result to {args.outfile}")
with open(args.outfile, "w") as f:
    f.write(result)

print("Done")
