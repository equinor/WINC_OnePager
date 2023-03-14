#!/usr/bin/env python
#GaP_17062022.py + SATNUM =2 is added for linear relperms inside the pipe (pipe_builder function) 
#GaP_17062022.py + well sketch visualization is added for QC 
#GaP_02052022.py  + the ACTNUM inside CARFIN is changed to TRANX/Y/Z for sake of runtime and visualization 
#GaP_05042022_FPT.py + the grids around legacy well inside LGR are de-activated 
#GaP.py + the reservoir is isolated from OVB 
# the only communication path within OVB and reservoir is pipe (ACTNUM inside CARFIN) 

from this import d
import numpy as np
from ecl.eclfile import EclFile
from ecl.eclfile import EclInitFile, EclRestartFile
from ecl.grid import EclGrid
from datetime import datetime  
from matplotlib.pyplot import show 
import matplotlib.pyplot as plt
#import rips 
import math
from os.path import exists

## making the 3D grid  

def make_stat_3d_main(vec,default_val = -1):

    ix = grid.export_index()
    aix = np.where(ix.active>=0)
    aix = np.where(init['PORV'][0].numpy_view()>0.0)
    nx,ny,nz = grid.get_dims()[:3]
    n= nx*ny*nz
    vec1d = default_val*np.ones(shape=(n,))
    vec1d[aix] = init[vec]
    vec3d = vec1d.reshape((nx,ny,nz), order='F')
    return vec3d

####################Loading the model############################## 
#simulation case without legacy well 

path = '/project/SCS/LEGACY_WELLS/ressim/Paris/PFT/'
simcase = 'GEN_NOLGR_T5'

grid = EclGrid(simcase + ".EGRID") 
#init = EclGrid(simcase + ".INIT") 
init = EclInitFile(grid ,simcase + ".INIT")
###############Main grid input data################### 

main_grd_i = 10
main_grd_j = 10
main_grd_min_k = 1
main_grd_max_k = 50

main_grd_dx =       make_stat_3d_main('DX',0)[main_grd_i-1,main_grd_j-1,main_grd_min_k-1]       #656.168 ft 
main_grd_dy =       make_stat_3d_main('DY',0)[main_grd_i-1,main_grd_j-1,main_grd_min_k-1]       #656.168 ft 
main_grd_dz_in_OB = make_stat_3d_main('DZ',0)[main_grd_i-1,main_grd_j-1,main_grd_min_k-1]  #500 ft 

no_of_layers_in_OB = 10   #Number of layer in OverBurden 
ref_depth = make_stat_3d_main('DEPTH',0)[main_grd_i-1,main_grd_j-1,main_grd_min_k-1] #Reference depth where LGR (legacy well starts) starts

#LRG name 
LGR_NAME = 'LEG_HIRES'



####################################Casing & Cement & Barrier input Data ########################

################################Conduction Casing Design  
cond_Casing_ID = 0.766  #30     inch casing 

cond_Casing_strt_depth = 40
cond_Casing_end_depth= 450

cond_oph_strt_depth = 4 
cond_oph_end_depth = 450

cond_Casing_strt_depth_cement = 180
cond_Casing_end_depth_cement= 450

###############################Surface Casing Design 
surf_Casing_ID = 0.339  #13-3/8 inch casing 

surf_Casing_strt_depth= 290
surf_Casing_end_depth= 726

surf_oph_strt_depth = 290
surf_oph_end_depth = 726

surf_Casing_strt_depth_cement= 726-(726-290)/3
surf_Casing_end_depth_cement= 726

###############################Production Casing Design 
prod_Casing_ID = 0.244  #9-5/8  inch casing 

prod_Casing_strt_depth = 650
prod_Casing_end_depth = 1450

prod_oph_strt_depth = 650
prod_oph_end_depth = 1600

prod_Casing_strt_depth_cement= 1450-(1450-500)/3
prod_Casing_end_depth_cement= 1450
######################################Barrier Desing 
has_barrier = 'yes'   # if the legacy well has barrier: 'yes' otherwise 'no'
barrier_ID = 0.244
barrier_strt_depth = 800
barrier_end_depth = 850 
barrier_perm = 0.01

#Permeability of the Tube and cement  
pipe_perm = 10000
cement_perm = 0.01 


if exists(LGR_NAME+ '.grdecl') == True:
    O = open(LGR_NAME+'.grdecl',"r+")
else: 
    O = open(LGR_NAME+'.grdecl',"x")

O.truncate(0)

#Calculating the size of grids for cement bond 
cond_bond = cond_Casing_ID - 2*surf_Casing_ID
surf_bond = surf_Casing_ID - prod_Casing_ID
prod_bond = surf_Casing_ID - prod_Casing_ID
#conversion of Casing diameter to cartesian system 
case_dim = [np.sqrt(0.25*np.pi*cond_Casing_ID**2),np.sqrt(0.25*np.pi*surf_Casing_ID**2) , np.sqrt(0.25*np.pi*prod_Casing_ID**2) ]
bond_dim = [cond_bond,surf_bond, prod_bond ] 
if round (min(case_dim + bond_dim),2) < 0.05 or round(min(case_dim + bond_dim),2) > 0.25:
    min_grd_size = 0.1 #round (min(case_dim + bond_dim),2)
else: 
    min_grd_size = round (min(case_dim + bond_dim),2)

#number of refined grids in tubes and cement bonds 
no_grd_surf_case = math.ceil (surf_Casing_ID/min_grd_size)
no_grd_cond_case = math.ceil (cond_Casing_ID/min_grd_size)
no_grd_prod_case = math.ceil (prod_Casing_ID/min_grd_size)

no_grd_surf_bond = math.ceil (surf_bond/min_grd_size)
no_grd_cond_bond = math.ceil (cond_bond/min_grd_size)
no_grd_prod_bond = math.ceil  (prod_bond/min_grd_size)

no_latral_grd_surf = no_grd_surf_case + no_grd_surf_bond*2

no_latral_grd_cond = no_grd_cond_case + no_grd_cond_bond*2

no_latral_grd_prod = no_grd_prod_case + no_grd_prod_bond*2

no_latral_fine_grd = max (no_latral_grd_prod,no_latral_grd_cond,no_latral_grd_surf)
LGR_size_fine_grd = [min_grd_size]*no_latral_fine_grd

LGR_size_no_outer =[min_grd_size*100] +[min_grd_size*10] +[min_grd_size]*no_latral_fine_grd+ [min_grd_size*10] +[min_grd_size*100]

LGR_sizes_xy = [(main_grd_dx- sum (LGR_size_no_outer))/2]+ LGR_size_no_outer + [(main_grd_dx- sum (LGR_size_no_outer))/2]

#all Z grids in OB is devided into 10 grids (hard coded), grids in reservoir remained unchanged 
#LGR_numb_z = (no_of_layers_in_OB-main_grd_min_k + 1) * [10] + (main_grd_max_k - main_grd_min_k + 1 -(no_of_layers_in_OB-main_grd_min_k+1) ) *[1] 
LGR_numb_z = (no_of_layers_in_OB-main_grd_min_k + 1) * [10] + (main_grd_max_k -no_of_layers_in_OB) *[1] 


LGR_numb_z = np.array (LGR_numb_z)

main_DZ = make_stat_3d_main('DZ',0)[main_grd_i-1,main_grd_j-1,(main_grd_min_k-1):main_grd_max_k]
LGR_size_z = np.divide(main_DZ,LGR_numb_z)
#ref_depth = make_stat_3d_main('DEPTH',0)[main_grd_i-1,main_grd_j-1,(main_grd_min_k-1)]
ref_depth = make_stat_3d_main('DEPTH',0)[main_grd_i-1,main_grd_j-1,(main_grd_min_k-1)] - (0.5*make_stat_3d_main('DZ',0)[main_grd_i-1,main_grd_j-1,(main_grd_min_k-1)])
LGR_intvl =  np.zeros ((LGR_numb_z.sum(),0))
LGR_depths= np.zeros ((LGR_numb_z.sum(),0))

LGR_index_z = np.arange(1,LGR_numb_z.sum()+1)



for i in range (0,LGR_numb_z.shape[0]): 
    LGR_intvl= np.append (LGR_intvl, np.repeat(LGR_size_z[i],LGR_numb_z[i]))  

   
# Depth conversion from field to metric. sample .EGRID model is in Field unit     
LGR_intvl [0] = LGR_intvl[0]+ ref_depth
LGR_depths = (np.cumsum(LGR_intvl))

main_DZ[0] = ref_depth + main_DZ[0]
MAIN_depths = (np.cumsum (main_DZ))


 # def pipe_builder (ID,strt_depth,end_depth,perm):
 #     if (LGR_depths[0] -strt_depth) > 0:
 #         k_min_pipe = 1 
 #     else:
 #        k_min_pipe = np.argmin(abs(LGR_depths-strt_depth))+ 1 

 #    if (LGR_depths[-1] -end_depth) < 0:
 #        k_max_pipe =  len(LGR_depths)
 #    else:
 #        k_max_pipe = np.argmin(abs(LGR_depths-end_depth))+ 1 
 #    #print (k_min_pipe, k_max_pipe)
    
 #    no_grd_pipe_x = math.floor (ID/min_grd_size)  
 #    x_min_pipe = math.ceil ((len (LGR_sizes_xy)-no_grd_pipe_x)/2 )
 #    if x_min_pipe < 0:
 #        print ('ERROR:The ID of tube is larger than refined area',file=O)
 #    x_max_pipe = x_min_pipe + no_grd_pipe_x 
 #    #print (x_min_pipe, x_max_pipe)
 
 #    no_grd_pipe_y = math.floor (ID/min_grd_size)  
 #    y_min_pipe = math.ceil((len (LGR_sizes_xy)-no_grd_pipe_x)/2)
 #    if y_min_pipe < 0:
 #        print ('ERROR:The ID of tube is larger than refined area',file=O)
 #    y_max_pipe = y_min_pipe + no_grd_pipe_y 
 #    #print (x_min_pipe, x_max_pipe) 
 #    print ('EQUALS',file=O)
 #    print ('--pipe with ID of',ID*39.37,'and perm of',perm,' were set in',LGR_NAME, 'Local Grid refinement',file=O)
 #    print ('PERMX','',perm,'',x_min_pipe,'',x_max_pipe,'',y_min_pipe,'',y_max_pipe,'',k_min_pipe,'',k_max_pipe,'','/',file=O)
 #    print ('PERMY','',perm,'',x_min_pipe,'',x_max_pipe,'',y_min_pipe,'',y_max_pipe,'',k_min_pipe,'',k_max_pipe,'','/',file=O)
 #    print ('PERMZ','',perm,'',x_min_pipe,'',x_max_pipe,'',y_min_pipe,'',y_max_pipe,'',k_min_pipe,'',k_max_pipe,'','/',file=O)
 #    print ('PORO','','0.99','',x_min_pipe,'',x_max_pipe,'',y_min_pipe,'',y_max_pipe,'',k_min_pipe,'',k_max_pipe,'','/',file=O)
 #    print ('--Transmisibilities of the edge of the pipe set to zero',file=O)
 #    print ('MULTX','',0,'',x_min_pipe,'',x_min_pipe,'',y_min_pipe,'',y_max_pipe,'',k_min_pipe,'',k_max_pipe,'','/',file=O)
 #    print ('MULTX','',0,'',x_max_pipe,'',x_max_pipe,'',y_min_pipe,'',y_max_pipe,'',k_min_pipe,'',k_max_pipe,'','/',file=O)
 #    print ('MULTX','',0,'',x_min_pipe,'',x_max_pipe,'',y_min_pipe,'',y_min_pipe,'',k_min_pipe,'',k_max_pipe,'','/',file=O)
 #    print ('MULTX','',0,'',x_min_pipe,'',x_max_pipe,'',y_max_pipe,'',y_max_pipe,'',k_min_pipe,'',k_max_pipe,'','/',file=O)

 #    print ('MULTY','',0,'',x_min_pipe,'',x_max_pipe,'',y_min_pipe,'',y_min_pipe,'',k_min_pipe,'',k_max_pipe,'','/',file=O)
 #    print ('MULTY','',0,'',x_min_pipe,'',x_max_pipe,'',y_max_pipe,'',y_max_pipe,'',k_min_pipe,'',k_max_pipe,'','/',file=O)
 #    print ('MULTY','',0,'',x_min_pipe,'',x_min_pipe,'',y_min_pipe,'',y_max_pipe,'',k_min_pipe,'',k_max_pipe,'','/',file=O)
 #    print ('MULTY','',0,'',x_max_pipe,'',x_max_pipe,'',y_min_pipe,'',y_max_pipe,'',k_min_pipe,'',k_max_pipe,'','/',file=O)

 #    print ('--New Saturation region for flow in pipe',file=O)
 #    print ('SATNUM','',2,'',x_min_pipe,'',x_max_pipe,'',y_min_pipe,'',y_max_pipe,'',k_min_pipe,'',k_max_pipe,'','/',file=O)
 #    print ('/',file=O)




def pipe_builder_with_oph_sec (ID,strt_depth,end_depth,openhole_strt_depth,openhole_end_depth,perm):
    if (LGR_depths[0] -strt_depth) > 0:
        k_min_pipe = 1 
    else:
        k_min_pipe = np.argmin(abs(LGR_depths-strt_depth))+ 1 

    if (LGR_depths[0] -openhole_strt_depth) > 0:
        k_min_hole = 1 
    else:
        k_min_hole = np.argmin(abs(LGR_depths-openhole_strt_depth))+ 1

    if (LGR_depths[-1] -end_depth) < 0:
        k_max_pipe =  len(LGR_depths)
    else:
        k_max_pipe = np.argmin(abs(LGR_depths-end_depth))+ 1 
    if (LGR_depths[-1] -openhole_end_depth) < 0:
        k_max_hole =  len(LGR_depths)
    else:
        k_max_hole = np.argmin(abs(LGR_depths-openhole_end_depth))+ 1 
    
    
    no_grd_pipe_x = math.floor ((ID-0.0001)/min_grd_size)  
    x_min_pipe = math.ceil ((len (LGR_sizes_xy)-no_grd_pipe_x)/2 )
    if x_min_pipe < 0:
        print ('ERROR:The ID of tube is larger than refined area',file=O)
    x_max_pipe = x_min_pipe + no_grd_pipe_x 
    #print (x_min_pipe, x_max_pipe)
 
    no_grd_pipe_y = math.floor ((ID-0.0001)/min_grd_size)  
    y_min_pipe = math.ceil((len (LGR_sizes_xy)-no_grd_pipe_x)/2)
    if y_min_pipe < 0:
        print ('ERROR:The ID of tube is larger than refined area',file=O)
    y_max_pipe = y_min_pipe + no_grd_pipe_y 
    #print (x_min_pipe, x_max_pipe) 
    print ('EQUALS',file=O)
    print ('--pipe with ID of',ID*39.37,'and perm of',perm,' were set in',LGR_NAME, 'Local Grid refinement',file=O)
    print ('PERMX','',perm,'',x_min_pipe,'',x_max_pipe,'',y_min_pipe,'',y_max_pipe,'',k_min_hole,'',k_max_hole,'','/',file=O)
    print ('PERMY','',perm,'',x_min_pipe,'',x_max_pipe,'',y_min_pipe,'',y_max_pipe,'',k_min_hole,'',k_max_hole,'','/',file=O)
    print ('PERMZ','',perm,'',x_min_pipe,'',x_max_pipe,'',y_min_pipe,'',y_max_pipe,'',k_min_hole,'',k_max_hole,'','/',file=O)
    print ('PORO','','0.99','',x_min_pipe,'',x_max_pipe,'',y_min_pipe,'',y_max_pipe,'',k_min_hole,'',k_max_hole,'','/',file=O)
    print ('SATNUM','',2,'',x_min_pipe,'',x_max_pipe,'',y_min_pipe,'',y_max_pipe,'',k_min_hole,'',k_max_hole,'','/',file=O)
    print ('--Transmisibilities of the edge of the pipe set to zero',file=O)
    print ('MULTX','',0,'',x_min_pipe,'',x_min_pipe,'',y_min_pipe,'',y_max_pipe,'',k_min_pipe,'',k_max_pipe,'','/',file=O)
    print ('MULTX','',0,'',x_max_pipe,'',x_max_pipe,'',y_min_pipe,'',y_max_pipe,'',k_min_pipe,'',k_max_pipe,'','/',file=O)
    print ('MULTX','',0,'',x_min_pipe,'',x_max_pipe,'',y_min_pipe,'',y_min_pipe,'',k_min_pipe,'',k_max_pipe,'','/',file=O)
    print ('MULTX','',0,'',x_min_pipe,'',x_max_pipe,'',y_max_pipe,'',y_max_pipe,'',k_min_pipe,'',k_max_pipe,'','/',file=O)

    print ('MULTY','',0,'',x_min_pipe,'',x_max_pipe,'',y_min_pipe,'',y_min_pipe,'',k_min_pipe,'',k_max_pipe,'','/',file=O)
    print ('MULTY','',0,'',x_min_pipe,'',x_max_pipe,'',y_max_pipe,'',y_max_pipe,'',k_min_pipe,'',k_max_pipe,'','/',file=O)
    print ('MULTY','',0,'',x_min_pipe,'',x_min_pipe,'',y_min_pipe,'',y_max_pipe,'',k_min_pipe,'',k_max_pipe,'','/',file=O)
    print ('MULTY','',0,'',x_max_pipe,'',x_max_pipe,'',y_min_pipe,'',y_max_pipe,'',k_min_pipe,'',k_max_pipe,'','/',file=O)

    print ('/',file=O)


def cement_bond_builder (ID, strt_depth, end_depth, perm):
    if (LGR_depths[0] -strt_depth) > 0:
        k_min_CB = 1 
    else:
        k_min_CB = np.argmin(abs(LGR_depths-strt_depth))+ 1 

    if (LGR_depths[-1] -end_depth) < 0:
        k_max_CB =  len(LGR_depths)
    else:
        k_max_CB = np.argmin(abs(LGR_depths-end_depth))+ 1 
    no_grd_pipe_x = math.floor (ID/min_grd_size)  
    x_min_pipe = math.ceil ((len (LGR_sizes_xy)-no_grd_pipe_x)/2)
    if x_min_pipe < 0:
        print ('ERROR:The ID of tube is larger than refined area',file=O)
    x_max_pipe = x_min_pipe + no_grd_pipe_x 
    #print (x_min_pipe, x_max_pipe)
 
    no_grd_pipe_y = math.floor (ID/min_grd_size)  
    y_min_pipe = math.ceil ((len (LGR_sizes_xy)-no_grd_pipe_x)/2)
    if y_min_pipe < 0:
        print ('ERROR:The ID of tube is larger than refined area',file=O)
    y_max_pipe = y_min_pipe + no_grd_pipe_y 
    print ('EQUALS',file=O)
    print ('--cement around pipe with ID of',ID*39.37,' and perm of',perm,' were set in',LGR_NAME,file=O)
    print ('--Top side',file=O)
    print ('PERMX','',perm,'',x_min_pipe-1,'',x_max_pipe+1,'',y_min_pipe-1,'',y_min_pipe-1,'',k_min_CB,'',k_max_CB,'','/',file=O)
    print ('PERMY','',perm,'',x_min_pipe-1,'',x_max_pipe+1,'',y_min_pipe-1,'',y_min_pipe-1,'',k_min_CB,'',k_max_CB,'','/',file=O)
    print ('PERMZ','',perm,'',x_min_pipe-1,'',x_max_pipe+1,'',y_min_pipe-1,'',y_min_pipe-1,'',k_min_CB,'',k_max_CB,'','/',file=O)
    print ('PORO','', 0.01,'',x_min_pipe-1,'',x_max_pipe+1,'',y_min_pipe-1,'',y_min_pipe-1,'',k_min_CB,'',k_max_CB,'','/',file=O)
    print ('',file=O)
    print ('--Bottom side',file=O)
    print ('PERMX','',perm,'',x_min_pipe-1,'',x_max_pipe+1,'',y_max_pipe+1,'',y_max_pipe+1,'',k_min_CB,'',k_max_CB,'','/',file=O)
    print ('PERMY','',perm,'',x_min_pipe-1,'',x_max_pipe+1,'',y_max_pipe+1,'',y_max_pipe+1,'',k_min_CB,'',k_max_CB,'','/',file=O)
    print ('PERMZ','',perm,'',x_min_pipe-1,'',x_max_pipe+1,'',y_max_pipe+1,'',y_max_pipe+1,'',k_min_CB,'',k_max_CB,'','/',file=O)
    print ('PORO','', 0.01,'',x_min_pipe-1,'',x_max_pipe+1,'',y_max_pipe+1,'',y_max_pipe+1,'',k_min_CB,'',k_max_CB,'','/',file=O)
    print ('',file=O)
    print ('--Left side',file=O) 
    print ('PERMX','',perm,'',x_min_pipe-1,'',x_min_pipe-1,'',y_min_pipe,'',y_max_pipe,'',k_min_CB,'',k_max_CB,'','/',file=O)
    print ('PERMY','',perm,'',x_min_pipe-1,'',x_min_pipe-1,'',y_min_pipe,'',y_max_pipe,'',k_min_CB,'',k_max_CB,'','/',file=O)
    print ('PERMZ','',perm,'',x_min_pipe-1,'',x_min_pipe-1,'',y_min_pipe,'',y_max_pipe,'',k_min_CB,'',k_max_CB,'','/',file=O)
    print ('PORO','', 0.01,'',x_min_pipe-1,'',x_min_pipe-1,'',y_min_pipe,'',y_max_pipe,'',k_min_CB,'',k_max_CB,'','/',file=O)
    print ('--Right side',file=O)
    print ('PERMX','',perm,'',x_max_pipe+1,'',x_max_pipe+1,'',y_min_pipe,'',y_max_pipe,'',k_min_CB,'',k_max_CB,'','/',file=O)
    print ('PERMY','',perm,'',x_max_pipe+1,'',x_max_pipe+1,'',y_min_pipe,'',y_max_pipe,'',k_min_CB,'',k_max_CB,'','/',file=O)
    print ('PERMZ','',perm,'',x_max_pipe+1,'',x_max_pipe+1,'',y_min_pipe,'',y_max_pipe,'',k_min_CB,'',k_max_CB,'','/',file=O)
    print ('PORO','', 0.01,'',x_max_pipe+1,'',x_max_pipe+1,'',y_min_pipe,'',y_max_pipe,'',k_min_CB,'',k_max_CB,'','/',file=O)
    print ('/',file=O)


def barrier_builder (ID,strt_depth, end_depth, perm):
    if (LGR_depths[0] -strt_depth) > 0:
        k_min_bar= 1 
    else:
        k_min_bar = np.argmin(abs(LGR_depths-strt_depth))+ 1 

    if (LGR_depths[-1] -end_depth) < 0:
        k_max_bar =  len(LGR_depths)
    else:
        k_max_bar = np.argmin(abs(LGR_depths-end_depth))+ 1 
    no_grd_bar_x = math.floor (ID/min_grd_size)  
    x_min_bar = math.ceil ((len (LGR_sizes_xy)-no_grd_bar_x)/2)
    if x_min_bar < 0:
        print ('ERROR:The ID of barrier is larger than refined area',file=O)
    x_max_bar = x_min_bar + no_grd_bar_x 
    
 
    no_grd_bar_y = math.floor (ID/min_grd_size)  
    y_min_bar = math.ceil ((len (LGR_sizes_xy)-no_grd_bar_x)/2)
    if y_min_bar < 0:
        print ('ERROR:The ID of tube is larger than refined area',file=O)
    y_max_bar = y_min_bar + no_grd_bar_y 
    print ('EQUALS',file=O)
    print ('--barrier with ID of',ID*39.37,'and perm of',perm,' were set in',LGR_NAME, file=O)
    print ('PERMX','',perm,'',x_min_bar,'',x_max_bar,'',y_min_bar,'',y_max_bar,'',k_min_bar,'',k_max_bar,'','/',file=O)
    print ('PERMY','',perm,'',x_min_bar,'',x_max_bar,'',y_min_bar,'',y_max_bar,'',k_min_bar,'',k_max_bar,'','/',file=O)
    print ('PERMZ','',perm,'',x_min_bar,'',x_max_bar,'',y_min_bar,'',y_max_bar,'',k_min_bar,'',k_max_bar,'','/',file=O)
    print ('PORO','','0.01','',x_min_bar,'',x_max_bar,'',y_min_bar,'',y_max_bar,'',k_min_bar,'',k_max_bar,'','/',file=O)
    print ('/',file=O)



#PRE-CARFIN, isolates OVB from reservoir 
print ('Prints isolating OVB from reservoir keywords in',LGR_NAME+'.grdecl file')   
print ('--isolating OVB from reservoir', file=O)
print ('EQUALS', file = O) 
print ('MULTX  0 ','1 ',grid.getNX(), '1 ',grid.getNY(),no_of_layers_in_OB+1,no_of_layers_in_OB+1,'/',file = O)
print ('MULTY  0 ','1 ',grid.getNX(), '1 ',grid.getNY(),no_of_layers_in_OB+1,no_of_layers_in_OB+1,'/',file = O)
print ('MULTZ  0 ','1 ',grid.getNX(), '1 ',grid.getNY(),no_of_layers_in_OB+1,no_of_layers_in_OB+1,'/',file = O)
print ('MULTX  1 ',main_grd_i, main_grd_i,main_grd_j, main_grd_j, no_of_layers_in_OB+1 ,no_of_layers_in_OB+1 ,'/',file =O)
print ('MULTY  1 ',main_grd_i, main_grd_i,main_grd_j, main_grd_j, no_of_layers_in_OB+1 ,no_of_layers_in_OB+1 ,'/',file =O)
print ('MULTZ  1 ',main_grd_i, main_grd_i,main_grd_j, main_grd_j, no_of_layers_in_OB+1 ,no_of_layers_in_OB+1 ,'/',file =O)
print ('/', file = O) 
print (' ', file = O) 
 

# prints CARFIN KEYWORD in Eclipse 
print ('Prints CARFIN Keywords in',LGR_NAME,'.grdecl file')
print ('CARFIN',file=O)
print (LGR_NAME ,main_grd_i, main_grd_i, main_grd_j, main_grd_j, main_grd_min_k, main_grd_max_k, len(LGR_sizes_xy), len(LGR_sizes_xy), sum(LGR_numb_z) , '/',file=O)

print (' ',file=O)
print ('NXFIN',file=O)
print (len(LGR_sizes_xy),'/',file=O)

print (' ',file=O)
print ('NYFIN',file=O)
print (len(LGR_sizes_xy),'/',file=O)

print (' ',file=O)
print ('NZFIN',file=O)
print (*LGR_numb_z, '/',file=O)
print (' ',file=O)
print ('HXFIN',file=O) 
print  (*[round (x/ min_grd_size,2) for x in LGR_sizes_xy],'/',file=O) 
print (' ',file=O) 
print ('HYFIN',file=O)
print  (*[round (x/ min_grd_size,2) for x in LGR_sizes_xy],'/',file=O) 
print (' ',file=O)
print (' ',file=O)

print ('Prints Casings, cement bonds and barrie(s) in',LGR_NAME+'.grdecl file') 
pipe_builder_with_oph_sec (cond_Casing_ID,cond_Casing_strt_depth,cond_Casing_end_depth,cond_oph_strt_depth,cond_oph_end_depth ,pipe_perm )
cement_bond_builder (cond_Casing_ID,cond_Casing_strt_depth_cement,cond_Casing_end_depth_cement,cement_perm )


pipe_builder_with_oph_sec (surf_Casing_ID,surf_Casing_strt_depth,surf_Casing_end_depth,surf_oph_strt_depth,surf_oph_end_depth ,pipe_perm )
cement_bond_builder (surf_Casing_ID,surf_Casing_strt_depth_cement,surf_Casing_end_depth_cement,cement_perm )


pipe_builder_with_oph_sec (prod_Casing_ID,prod_Casing_strt_depth,prod_Casing_end_depth,prod_oph_strt_depth,prod_oph_end_depth ,pipe_perm )
cement_bond_builder (prod_Casing_ID,prod_Casing_strt_depth_cement,prod_Casing_end_depth_cement,cement_perm )

##optional: to add barrier inside the well
if has_barrier == 'yes':
    barrier_builder (barrier_ID, barrier_strt_depth,barrier_end_depth,barrier_perm) 

print ('Prints isolating OVB from reservoir in the LGR in ',LGR_NAME+'.grdecl file') 
#Trans. modification to isolate OVB from Reservoir inside the LGR
print ('--isolating OVB from reservoir in the LGR',file = O)
print ('EQUALS',file = O ) 
print ('MULTX  0 ','1 ',len(LGR_sizes_xy), '1 ',len(LGR_sizes_xy),(no_of_layers_in_OB-main_grd_min_k+1)*10+1,(no_of_layers_in_OB-main_grd_min_k+1)*10+1,'/',file = O )
print ('MULTY  0 ','1 ',len(LGR_sizes_xy), '1 ',len(LGR_sizes_xy),(no_of_layers_in_OB-main_grd_min_k+1)*10+1,(no_of_layers_in_OB-main_grd_min_k+1)*10+1,'/',file = O )
print ('MULTZ  0 ','1 ',len(LGR_sizes_xy), '1 ',len(LGR_sizes_xy),(no_of_layers_in_OB-main_grd_min_k+1)*10+1,(no_of_layers_in_OB-main_grd_min_k+1)*10+1,'/',file = O )


reopen_ID = min ([cond_Casing_ID, surf_Casing_ID, prod_Casing_ID])
no_grd_reopen_x = math.floor (reopen_ID/min_grd_size)  
x_min_reopen = math.ceil ((len (LGR_sizes_xy)-no_grd_reopen_x)/2)
x_max_reopen = x_min_reopen + no_grd_reopen_x 

no_grd_reopen_y = math.floor (reopen_ID/min_grd_size)  
y_min_reopen = math.ceil ((len (LGR_sizes_xy)-no_grd_reopen_y)/2)
y_max_reopen = y_min_reopen + no_grd_reopen_y 

print ('MULTZ 1 ', x_min_reopen-1,x_max_reopen+1,y_min_reopen-1 ,y_max_reopen+1 ,(no_of_layers_in_OB-main_grd_min_k+1)*10+1 ,(no_of_layers_in_OB-main_grd_min_k+1)*10+1 ,'/',file = O)

print ('/', file = O) 

print ('ENDFIN', file = O) 
O.close()

#########################Well Sketch Design########################3
print ('Visualizes the well and its locationin OVB and reservoir...') 

def pipe_plotter (ID, casing_strt_depth, casing_end_depth,casing_strt_depth_cement, casing_end_depth_cement, oph_strt_depth, oph_end_depth):
    #Casing lines 
    plt.vlines (x = -ID/2 , ymin = - casing_end_depth , ymax = -casing_strt_depth ) 
    plt.vlines (x = (-ID/2)+ID , ymin = -casing_end_depth , ymax = -casing_strt_depth )
    plt.scatter(-ID/2,-casing_end_depth, marker=8,color = 'black')
    plt.scatter(ID/2,-casing_end_depth, marker=9, color = 'black')
    if (LGR_depths[0] -casing_strt_depth) > 0:
        k_min_pipe = 1 
    else:
        k_min_pipe = np.argmin(abs(LGR_depths-casing_strt_depth))+ 1 

    if (LGR_depths[0] -oph_strt_depth) > 0:
        k_min_hole = 1 
    else:
        k_min_hole = np.argmin(abs(LGR_depths-oph_strt_depth))+ 1

    if (LGR_depths[-1] -casing_end_depth) < 0:
        k_max_pipe =  len(LGR_depths)
    else:
        k_max_pipe = np.argmin(abs(LGR_depths-casing_end_depth))+ 1 
    if (LGR_depths[-1] -oph_end_depth) < 0:
        k_max_hole =  len(LGR_depths)
    else:
        k_max_hole = np.argmin(abs(LGR_depths-oph_end_depth))+ 1 
    
    plt.annotate( str(k_max_hole), (-ID/2-0.5,-oph_end_depth))
    plt.annotate( str(k_min_hole), (-ID/2-0.5,-oph_strt_depth))

    plt.annotate( str(k_max_pipe), (ID/2+0.5,-casing_end_depth))
    plt.annotate( str(k_min_pipe), (ID/2+0.5,-casing_strt_depth))
    plt.annotate( 'Openhole K index', (-3.5,5))
    plt.annotate( 'Casing K index', (+2,0))
    #cement bonds 
    plt.fill_between ([-ID/2, -ID/2-0.08], -casing_end_depth_cement,- casing_strt_depth_cement, color='grey', alpha = 0.5  )  #cement bond thickness is hardcoded to 0.05 for visualization 
    plt.fill_between ([ID/2,ID/2+0.08 ], -casing_end_depth_cement,-casing_strt_depth_cement, color='grey', alpha = 0.5  )        #cement bond thickness is hardcoded to 0.05 for visualization 
    #fill in pipe 
    plt.fill_between(  [-ID/2,(-ID/2)+ID] , -oph_end_depth ,- oph_strt_depth, color = 'green', alpha = 0.2 )

def barrier_plotter (ID, barrier_strt_depth,  barrier_end_depth):
    plt.fill_between(  [-ID/2,(-ID/2)+ID] , - barrier_end_depth ,- barrier_strt_depth, color = 'red', alpha = 0.5 )




pipe_plotter (cond_Casing_ID, cond_Casing_strt_depth, cond_Casing_end_depth, cond_Casing_strt_depth_cement, cond_Casing_end_depth_cement,  cond_oph_strt_depth, cond_oph_end_depth)
pipe_plotter (surf_Casing_ID, surf_Casing_strt_depth, surf_Casing_end_depth, surf_Casing_strt_depth_cement, surf_Casing_end_depth_cement,  surf_oph_strt_depth, surf_oph_end_depth)
pipe_plotter (prod_Casing_ID, prod_Casing_strt_depth, prod_Casing_end_depth, prod_Casing_strt_depth_cement, prod_Casing_end_depth_cement,  prod_oph_strt_depth, prod_oph_end_depth)
barrier_plotter(barrier_ID,  barrier_strt_depth, barrier_end_depth )  


#Overburden/Reservoir Design 
plt.fill_between(  [-4,4] , -make_stat_3d_main('DEPTH',0)[main_grd_i-1,main_grd_j-1,0], - make_stat_3d_main('DEPTH',0)[main_grd_i-1,main_grd_j-1,no_of_layers_in_OB] , color = 'brown', alpha = 0.1 )
plt.fill_between(  [-4,4] , -make_stat_3d_main('DEPTH',0)[main_grd_i-1,main_grd_j-1,no_of_layers_in_OB], - make_stat_3d_main('DEPTH',0)[main_grd_i-1,main_grd_j-1,59] , color = 'brown', alpha = 0.5 )
#plt.plot(np.log10(make_stat_3d_main('PERMZ',0)[main_grd_i-1,main_grd_j-1,0:139]),-make_stat_3d_main('DEPTH',0)[main_grd_i-1,main_grd_j-1,0:139])
plt.xlim(-4,4)
plt.show()  

