""" generate grdecl grid
"""
from typing import List, Union

import os
import math
from typing import TextIO

import numpy as np

from .CARFIN_builders import (
    CARFIN_pipe_and_cement_bond_builder,
    CARFIN_oph_builder,
    CARFIN_barrier_builder,
)

from ..models import (
    PipeCementModel,
    ElemModel,
)

#
def pre_CARFIN(LGR_NAME: str,
               NX: int, NY: int,
               main_grd_i: int, main_grd_j: int, 
               no_of_layers_in_OB: int, 
               O: TextIO):
    """ preCARFIN

        Before introducing the CARFIN, we need to make sure there is a sealing layer between the overburden and the reservoir.
        
        We borrow  the topmost layer of the reservoir and make it seal against both vertical (MULTZ=0) and side (MULTX, MULTY= 0) 
        flow toward the overburden. 

        But, we still need to make sure that the topmost layer is open in where the LGR is located to allow the only flow communication 
        between the overburden and the reservoir being through the well 

        Later we will close the rest of topmost layer inside the LGR and only open an area of LGR where the main pipe is paths.

        Depend on the case, we may want to open the cement bonds in that area.   

        Args:

            LGR_NAME (str): output filename, without the suffix '.grdecl'
            NX (int): number of grids in x direction (coarse grid)
            NY (int): number of grids in y direction (coarse grid)
            main_grd_i (int): x index of well location in coarse grid
            main_grd_j (int): y index of well location in coarse grid
            no_of_layers_in_OB (int): number of layers in overburden (coarse grid)
            O (TextIO): opened file handle            
    """

    # PRE-CARFIN, isolates OVB from reservoir 
    print ('...Prints isolating OVB from reservoir keywords in', LGR_NAME+'.grdecl file')

    print ('--isolating OVB from reservoir', file=O)
    print ('EQUALS', file = O) 
    print ('MULTX  0 ','1 ',NX, '1 ',NY,no_of_layers_in_OB+1,no_of_layers_in_OB+1,'/',file = O)
    print ('MULTY  0 ','1 ',NX, '1 ',NY,no_of_layers_in_OB+1,no_of_layers_in_OB+1,'/',file = O)
    print ('MULTZ  0 ','1 ',NX, '1 ',NY,no_of_layers_in_OB+1,no_of_layers_in_OB+1,'/',file = O)
    print ('MULTX  1 ',main_grd_i, main_grd_i,main_grd_j, main_grd_j, no_of_layers_in_OB+1, no_of_layers_in_OB+1, '/',file =O)
    print ('MULTY  1 ',main_grd_i, main_grd_i,main_grd_j, main_grd_j, no_of_layers_in_OB+1, no_of_layers_in_OB+1, '/',file =O)
    print ('MULTZ  1 ',main_grd_i, main_grd_i,main_grd_j, main_grd_j, no_of_layers_in_OB+1, no_of_layers_in_OB+1, '/',file =O)
    print ('/', file = O) 
    print (' ', file = O) 
    
def CARFIN_keywords(LGR_NAME: str,
                    main_grd_i: int, main_grd_j: int, 
                    main_grd_min_k: int, main_grd_max_k: int, 
                    LGR_sizes_xy: List[float], 
                    LGR_numb_z: np.ndarray, 
                    min_grd_size: float,
                    O: TextIO):
    """ CARFIN main parameters

        CARFIN is the main keyword which introduces LGR into the  grid. 

        The keyword has the following parameters which had be filled in: 

            LGR Name: the name of LGR that we want to create. 
                      This name is the same as the name that the script uses to generate the output file, i.e. LGRNAME.grdecl 

            There are 6 parameter which define where the LGR should be located in the main grid which are:

                main grid i-start, main grid i-end, main grid j-start, main grid j-end, main grid k-start, main grid k-end. 

            In our case, we always have a vertical well, therefore main grid i-start = main grid i-end are the same for j direction.  
    
            NXFIN and NYFIN: they dictate how many times per grid (in the same order) should  be chopped. In our case, 
                             there is only one grid assigned in the X and Y directions. Therefore, the numbers will be the same as 
                             8th and 9th parameters in CARFIN. 

            NZFIN: unlike X,Y direcitons, there are multiple grids involved in the Z direction of the main grid. 
                    The main grids in the overburden in the Z direction are hard-coded to be chopped to 10 grids.
                    The main grids in the Z direction in the reservoir should remain un-refined. Therefore, the assigned number of chops should be equal to 1. 
                    This is all already stored in the LGR_numb_z.

            HXFIN and HYFIN: They are going to provide multiple sizes in X and Y direction. This is mainly high resolution down to min_grd_size 
                             in the center of the grid and increasing the size of the grids toward the outermost side of the grid. 

                            The HXFIN and HYFIN are working not with real sizes but with “ratio”. 
                            The lowest grid size has the ratio equal to min_grd_size/min_grd_size = 1. 

                            Then the rest of grid sizes are calculated as follows:
                                Grid size of the refined grid /min_grd_size. 
                                The min_grid_size presented in slide7   
                                Grid sizes in x,y direction are presented in the slide 8

        Args:

            LGR_NAME (str): output filename, without the suffix '.grdecl'
            main_grd_i (int): x index of well location in coarse grid
            main_grd_j (int): y index of well location in coarse grid
            main_grd_min_k (int): minimum z index in coarse grid
            main_grd_max_k (int): maximum z index in coarse grid
            LGR_sizes_xy (list[float]): LGR xy grid intervals
            LGR_numb_z (np.ndarray): number of chops for each main DZ
            min_grd_size (float): minimum grid size
            O (TextIO): opened file handle            
    """

    # prints CARFIN KEYWORD in Eclipse 
    print ('...Prints CARFIN Keywords in', LGR_NAME+'.grdecl file')

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

def coreCARFIN(LGR_NAME: str,
                casing_list: List[Union[PipeCementModel, ElemModel]],
                barrier_list: List[ElemModel],
                LGR_sizes_xy: List[float], 
                LGR_depths: np.ndarray, 
                min_grd_size: float,
                O: TextIO):
    """ CARFIN for the main elements

        Args:

            LGR_NAME (str): output filename, without the suffix '.grdecl'
            casing_list (list[PipeCementModel]): contains list of casing geometry
            barrier_list (list[ElemModel]): contains list of barriers
            LGR_sizes_xy (list[float]): LGR xy grid intervals
            LGR_depths (np.ndarray): specific size of each LGR grid in the Z direction
            min_grd_size (float): minimum grid size
            O (TextIO): opened file handle            

    """

    print ('...Prints Casings, cement bonds and barrie(s) in', LGR_NAME+'.grdecl file')

    # loop around Conductor casing, Surface casing, and Production casing
    for idx, casing in enumerate(casing_list):
        if type(casing) is PipeCementModel:
            CARFIN_pipe_and_cement_bond_builder(casing, 
                                                LGR_sizes_xy, LGR_depths, 
                                                min_grd_size,
                                                LGR_NAME, 
                                                O)
        elif type(casing) is ElemModel:
            CARFIN_oph_builder(casing, 
                                LGR_sizes_xy, LGR_depths, 
                                min_grd_size,
                                LGR_NAME, 
                                O)
        else:
            continue
            
    # 4. optional: to add barrier inside the well
    if len(barrier_list):

        for _, barrier in enumerate(barrier_list):

            CARFIN_barrier_builder(barrier, 
                                    LGR_sizes_xy, LGR_depths, 
                                    min_grd_size,
                                    LGR_NAME, 
                                    O)

def endCARFIN(LGR_NAME: str,
             reopen_ID: float,
             LGR_sizes_xy: List[float], 
             main_grd_min_k: float, 
             min_grd_size: float,
             no_of_layers_in_OB: int,
             O):
    """ ENDFIN

        Before we finish the CARFIN with ENDFIN, we will have to close the area around the well inside the LGR:

            As mentioned in the slide 22, we opened the whole LGR in the topmost layer of the reservoir (and closed the rest of topmost layer.

            In this stage we will have to make sure there is no flow around the well inside the LGR from the reservoir to the overburden. 

            To do that we will first close the whole topmost layer inside the LGR:

                The whole LGR grids in the x,y directions should be closed in the topmost layer but inside LGR. 

                Then we will  set the MULTX, MULTY and MULTZ multipliers to 0 to make sure the area is whole closed.

            Then we re-open the area where the well paths. In order to do that, we will have to (again) find that which pipe (with which ID) 
            is pathing through the caprock (from overburden to the caprock):

                The code always assumes the narrowest pipe passing through the caprock.  

                We find X_min, X_max, Y_min, Y_max of the narrowest pipe.
                
            Knowing, how many layer are allocated to overburden and given 10 chops per grid in the Z direction, we can find the k_min and k_max of 
            the pipe where we want to leave it open. 
  
            Then, we let the MULTZ (only) equal to one in that domain. We don't want to allow side movement of CO2 in the caprock. 

        Args:

            LGR_NAME (str): output filename, without the suffix '.grdecl'
            reopen_ID (str): minimum ID of all casings
            LGR_sizes_xy (list[float]): LGR xy grid intervals
            main_grd_min_k (int): minimum z index in coarse grid
            no_of_layers_in_OB (int): number of layers in overburden (coarse grid)
            O (TextIO): opened file handle                 
    """

    # close the whole topmost layer inside the LGR
    print ('...Prints isolating OVB from reservoir in the LGR in', LGR_NAME+'.grdecl file')

    #Trans. modification to isolate OVB from Reservoir inside the LGR
    print ('--isolating OVB from reservoir in the LGR',file = O)
    print ('EQUALS',file = O ) 
    print ('MULTX  0 ','1 ',len(LGR_sizes_xy), '1 ',len(LGR_sizes_xy),(no_of_layers_in_OB-main_grd_min_k+1)*10+1,(no_of_layers_in_OB-main_grd_min_k+1)*10+1,'/',file = O )
    print ('MULTY  0 ','1 ',len(LGR_sizes_xy), '1 ',len(LGR_sizes_xy),(no_of_layers_in_OB-main_grd_min_k+1)*10+1,(no_of_layers_in_OB-main_grd_min_k+1)*10+1,'/',file = O )
    print ('MULTZ  0 ','1 ',len(LGR_sizes_xy), '1 ',len(LGR_sizes_xy),(no_of_layers_in_OB-main_grd_min_k+1)*10+1,(no_of_layers_in_OB-main_grd_min_k+1)*10+1,'/',file = O )

    # re-open the area where the well paths

    # x
    no_grd_reopen_x = math.floor(reopen_ID/min_grd_size)
    x_min_reopen = math.ceil((len (LGR_sizes_xy) - no_grd_reopen_x)/2)
    x_max_reopen = x_min_reopen + no_grd_reopen_x

    # y
    no_grd_reopen_y = math.floor(reopen_ID/min_grd_size)  
    y_min_reopen = math.ceil((len (LGR_sizes_xy) - no_grd_reopen_y)/2)
    y_max_reopen = y_min_reopen + no_grd_reopen_y 

    print ('MULTZ 1 ', x_min_reopen-1,x_max_reopen+1,y_min_reopen-1 ,y_max_reopen+1 ,(no_of_layers_in_OB-main_grd_min_k+1)*10+1 ,(no_of_layers_in_OB-main_grd_min_k+1)*10+1 ,'/',file = O)

    print ('/', file = O) 

    print ('ENDFIN', file = O)
