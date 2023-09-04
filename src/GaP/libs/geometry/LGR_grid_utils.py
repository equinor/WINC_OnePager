
import math
import numpy as np

def compute_cement_bonds(cond_Casing_ID: float, 
                         surf_Casing_ID: float, 
                         prod_Casing_ID: float) -> list[float, float, float]:
    """ compute cement bonds
    """

    # Calculating the size of grids for cement bond 
    # TODO(hzh): why multiply 2?
    cond_bond = cond_Casing_ID - surf_Casing_ID*2
    surf_bond = surf_Casing_ID - prod_Casing_ID
    prod_bond = surf_Casing_ID - prod_Casing_ID

    return [cond_bond, surf_bond, prod_bond]

def compute_min_grid_size(cond_Casing_ID: float, 
                          surf_Casing_ID: float, 
                          prod_Casing_ID: float) -> float:
    """ compute minimum grid size 

    The min. sizes happen when pipes are going into each other:

        From D&W, we know some cosignings (if not all) have overlap with each other and that's why we have annules flow. 

        In addition, the thickness of the cement bond could be the difference between the casings. 

        Finally, we don't want to go below 5cm grids, and above 25 cm is unphysically too large. Then we hard-code 
        the minim grid size to 10 cm in those case,
        
        Otherwise, let the code decide on the minimum grid size that we want to have. 

    Args:
        cond_Casing_ID (float): ID of conductor casing
        surf_Casing_ID (float): ID of surface casing
        prod_Casing_ID (float): ID of production casing

    Returns:
        float: minimum grid size
    """

    # Conversion of Casing diameter to cartesian system 
    # dimensions for casings
    case_dim = [np.sqrt(0.25*np.pi*cond_Casing_ID**2),
                np.sqrt(0.25*np.pi*surf_Casing_ID**2), 
                np.sqrt(0.25*np.pi*prod_Casing_ID**2)]


    # dimensions for cement bonds
    bond_dim = compute_cement_bonds(cond_Casing_ID, 
                                    surf_Casing_ID, 
                                    prod_Casing_ID)

    # compute the min grid size
    if round (min(case_dim + bond_dim), 2) < 0.05 or round(min(case_dim + bond_dim), 2) > 0.25:
        min_grd_size = 0.1 #round (min(case_dim + bond_dim),2)
    else: 
        min_grd_size = round (min(case_dim + bond_dim), 2)

    return min_grd_size

def compute_max_num_of_fine_grid_xy(cond_Casing_ID: float, 
                                    surf_Casing_ID: float, 
                                    prod_Casing_ID: float, 
                                    min_grd_size: float) -> int:
    """ compute maximum number of the fine grid in x-y direction

        For each casing we have its ID , therefore, the number of grids for each casing we should have round up casing_ID/min_grd_size.

        From the thickness of the cement bond, accordingly we can assume how many grids do we want that for this. 

        Then, per casing the total number of grids = number of grids to mimic the pipe + 2 x number of grids for the cement bond. 

        Args:
            cond_Casing_ID (float): ID of conductor casing
            surf_Casing_ID (float): ID of surface casing
            prod_Casing_ID (float): ID of production casing
            min_grd_size (float): minimum grid size

        Returns:
            LGR_sizes_xy (list[float]): LGR xy grid intervals
    """

    # Calculating the size of grids for cement bond 
    cond_bond, surf_bond, prod_bond = compute_cement_bonds(cond_Casing_ID, 
                                                            surf_Casing_ID, 
                                                            prod_Casing_ID)

    # Number of refined grids in tubes and cement bonds

    # 1. tubes
    no_grd_cond_case = math.ceil(cond_Casing_ID/min_grd_size)
    no_grd_surf_case = math.ceil(surf_Casing_ID/min_grd_size)
    no_grd_prod_case = math.ceil(prod_Casing_ID/min_grd_size)

    # 2. cement bonds
    no_grd_cond_bond = math.ceil(cond_bond/min_grd_size)
    no_grd_surf_bond = math.ceil(surf_bond/min_grd_size)
    no_grd_prod_bond = math.ceil(prod_bond/min_grd_size)

    # 3. TODO(hzh): why add them together and why multiply 2?
    no_latral_grd_cond = no_grd_cond_case + no_grd_cond_bond*2
    no_latral_grd_surf = no_grd_surf_case + no_grd_surf_bond*2
    no_latral_grd_prod = no_grd_prod_case + no_grd_prod_bond*2

    # print('ooooo', cond_bond, surf_bond, prod_bond)
    # print('xxxxx', no_latral_grd_cond, no_latral_grd_surf, no_latral_grd_prod)
    # print('yyyyy', no_grd_cond_case, no_grd_surf_case, no_grd_prod_bond)

    # 4. find maximum of the three
    no_latral_fine_grd = max(no_latral_grd_cond, no_latral_grd_surf, no_latral_grd_prod)

    return no_latral_fine_grd

def compute_LGR_xy(no_latral_fine_grd: int, 
                   main_grd_dx: float, 
                   min_grd_size: float) -> list[float]:
    """ LGR grid list in x-y direction

        Thinking laterally we need ones to the right and one to the left

        Args:
            no_latral_fine_grd (int): number of LGR lateral fine grids, including casing and cement bond
            main_grd_dx (float): DX value on coarse grid at the well location
            min_grd_size (float): minimum grid size

        Returns:
            LGR_sizes_xy (list[float]): LGR xy grid interval values
    """

    # # the fine grid
    # LGR_size_fine_grd = [min_grd_size]*no_latral_fine_grd

    # 5. the list
    # add transistion zone between coarse and fine grids
    LGR_size_no_outer = [min_grd_size*100] + \
                        [min_grd_size*10] + \
                        [min_grd_size] * no_latral_fine_grd + \
                        [min_grd_size*10] + \
                        [min_grd_size*100]

    # 6. LGR_sizes_xy, adding dx on both sides to make it one cell size
    LGR_sizes_xy = [(main_grd_dx - sum(LGR_size_no_outer))/2] + LGR_size_no_outer + [(main_grd_dx - sum(LGR_size_no_outer))/2]

    return LGR_sizes_xy

def compute_LGR_z(main_DZ,
                    ref_depth,
                    main_grd_min_k, main_grd_max_k, 
                    no_of_layers_in_OB) -> tuple[np.ndarray, np.ndarray]:
    """ LGR grid in z direction

    The idea is not to refine too much in the Z direction, mainly because of the computational time. 

    However, it we don't refine at all, the resolution in the Z direction in overburden might be too low. Then, it would limit 
    our access to a specific depth, when we want to set start and end depth of each casing. 

    First we define how many grids do we want? 

        We hard-coded the number of grids in overburden to 10 per layer. i.e. each layer is refined into 10 layers. 
        For example, if we have 5 coarse grids in the overburden, in the LGR we will have 5x10=50 layers. 
    
        Then we convert the list to numpy array to be able to do calculation on that.
    
        Later we will see the LGR_numb_z can directly be used in the NZFIN keyword under CARFIN  

    Args:
        main_DZ (np.ndarray): thickness or interval of each layer (coarse grid)
        main_grd_min_k (int): minimum z index (coarse grid)
        main_grd_max_k (int): maximum z index (coarse grid)
        no_of_layers_in_OB (int): number of layers in Overburden (coarse grid)

    Returns:

        LGR_depths (np.ndarray): depth values of each LGR z 
        LGR_numb_z (np.ndarray): number of refined layers in each coarse layer        
    """

    #all Z grids in OB is devided into 10 grids (hard coded), grids in reservoir remained unchanged

    #LGR_numb_z = (no_of_layers_in_OB-main_grd_min_k + 1) * [10] + (main_grd_max_k - main_grd_min_k + 1 -(no_of_layers_in_OB-main_grd_min_k+1) ) *[1] 

    # TODO(hzh): the number 10 is hard-coded?
    # refine only coarse grid in ovb, and keep reservoir grid unchanged
    LGR_numb_z = (no_of_layers_in_OB - main_grd_min_k + 1) * [10] + (main_grd_max_k - no_of_layers_in_OB) * [1] 
    LGR_numb_z = np.array(LGR_numb_z)

    #
    assert LGR_numb_z.size == main_DZ.size, 'dimension should match'
    # **** LGR_size_z
    LGR_size_z = np.divide(main_DZ, LGR_numb_z)

    # 
    LGR_intvl  = np.zeros((LGR_numb_z.sum(),0))
    # LGR_depths = np.zeros((LGR_numb_z.sum(),0))

    # LGR_index_z = np.arange(1, LGR_numb_z.sum()+1)

    for i in range (0, LGR_numb_z.shape[0]): 
        LGR_intvl= np.append(LGR_intvl, np.repeat(LGR_size_z[i], LGR_numb_z[i]))

    # Depth conversion from field to metric. sample .EGRID model is in Field unit
    # TODO(hzh): why add ref_depth?
    LGR_intvl [0] = LGR_intvl[0] + ref_depth
    LGR_depths = (np.cumsum(LGR_intvl))

    # main_DZ[0] = ref_depth + main_DZ[0]
    # MAIN_depths = (np.cumsum (main_DZ))

    return LGR_depths, LGR_numb_z
