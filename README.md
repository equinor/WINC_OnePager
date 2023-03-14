# GaP
A tool to build legacy well representations on reservoir simulation grids

<details><summary><h2>Introduction</h2></summary>

To be able to model a legacy well in reservoir scale, we need to make sure all of the elements including  multiple casings with different OD, cement bonds, barriers, the open area between barriers inside the casings, etc. are considered.
Then, the existing simple well models available in the commercial simulators (Eclipse, PFT, IX) are not able to include those details. They just introduce a node where the flow will be discharged from/to the grid to the node.

The way around it is to define the well as a part of the reservoir by manipulating the local grids (LGR) and properties of the grid, in that setting:

- **The pipes** can be mimicked by a very narrow (in the size of ID of a pipe-20 to 50 cm) and high permeability grids. Those grids could be surrounded by grids with zero transmissibility.
- **The open hole sections** can be modelled the same as pipe, but without zero transmissibility around the high perm. area to allow moving of fluids side ways from the piper. This is particularly important in the cases where there is a drilled but uncased hole in the setting.
- **The barriers** can be mimicked by very low permeability with the same size of pipe ID and have their own start and end depth.
- **The cement bond** can be mimicked by low perm vertical layer adjustment to the casing with given depth and interval.

Each unit above should have specified start and end depths.

In this context, the GaP script generates the LGR, changes the properties of the LGR to reflect all requested elements, and isolates reservoir from the overburden, opens the communication between reservroir and overburden only through the well. 
</details>

<details><summary><h2>GaP code structure</h2></summary>

The script is structured in the following main parts:

1. Input data
   1. The main grid
   2. Start, end and ID of each element.
   3. Permeability of elements (pipe_perm, barrier_perm)
2. LGR main part generation
   1. This tries to make a refine-enough grids in the middle of LGR, then the size of outer layers is increased by a factor of 10.
3. Functions which are element builders:
   1. Independent, explicit functions which can build pipe, cement bond and barriers.
   2. Those functions can be called several times, depend on how many pipes, or barriers or cement bonds do we have in the model. 
</details>
   
<details><summary><h2>Input data: <code>Main grid input data</code></h2></summary>

- In order to know the main-grid size, depths, etc. we need to load the main grid where there is no LGR is allocated  yet.
- The `make_stat_3d_main` function, is a standard way of loading a grid in 1D, exclude the in-active cells and change it to the 3D vector. Then i,j,k from the main grid can be called whenever needed and the function reports back the requested properties.

    ```python
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
    ```

- In all of the cases, LGR lies in one i,j and has start and end K values-the well supposed to be vertical.
    ```python
    main_grd_i = 10
    main_grd_j = 10
    main_grd_min_k = 1
    main_grd_max_k = 50
    ```
    
- To make it more automatic, the DX, DY,DZ (to be chopped into many grids) are read from the main grid.
    ```python
    main_grd_dx =       make_stat_3d_main('DX',0)[main_grd_i-1,main_grd_j-1,main_grd_min_k-1]
    main_grd_dy =       make_stat_3d_main('DY',0)[main_grd_i-1,main_grd_j-1,main_grd_min_k-1] 
    main_grd_dz_in_OB = make_stat_3d_main('DZ',0)[main_grd_i-1,main_grd_j-1,main_grd_min_k-1] 
    ```

- `No_of_layers_in_overburden` is important to know where the zero.-trans between reservoir and overburden should be located.
</details>

<details><summary><h2>Input data: <code>Casing & Cement & Barrier input Data</code></h2></summary>

In this section the specification of each casing is specified.
- The parameters ending with `_oph_strt_depth` and `_oph_end_depth` define the top and bottom of sections drilled but not necessarily cased. Similarly, for sections with casing, the `_Casing_strt_depth` and `_Casing_end_depth` are defined. Parameters ending with `_cement`define the top and bottom of the cement bond.
- The `_oph_strt_depth` should always be shallower or equal to `_Casing_start_depth`. And accordingly, the `_oph_end_depth` should always be deeper (or equal) than the `_Casing_end_depth`.



    ```python
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
    ```
    
</details>


<details><summary><h2>Output file generation</h2></summary>

- After finishing with the input parameters, an output file is created. 
- The output file outside the script environment is called the name of the CARFIN or LGR and ends with .grdecl 

For example: LGRNAME.grdecl 
- The output file internally is called “O”. 
- Later we will see that all of the prints goes directly into the “O”. 
- If the file is name is not generated before, the Code creates one for the user, 
- If the file is already generated by the previous runs, it will be cleaned (O.truncate(0) ) and will be re-written.
  
  ```python
  if exists(LGR_NAME+ '.grdecl') == True:
    O = open(LGR_NAME+'.grdecl',"r+")
  else: 
    O = open(LGR_NAME+'.grdecl',"x")

  O.truncate(0)
  ```
</details>


<details><summary><h2>Minimum grid size</h2></summary>

One of the most important parameters that we should find out from the beginning, is how refined our model should be. 
We should be looking for the minimin sizes that we are going to model. 
- The min. sizes happen when pipes are going into each other.
- From D&W, we know some cosignings (if not all) have overlap with each other and that’s why we have annules flow. 
- In addition, the thickness of the cement bond could be the difference between the casings. 
- Finally, we don’t want to go below 5cm grids, and above 25 cm is unphysically too large. Then we hard-code the minim grid size to 10 cm in those case, 
- Otherwise, let the code decide on the minimum grid size that we want to have. 
  
  ```python
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
  ```
</details>


<details><summary><h2>Number of fine grids in the whole LGR in X-Y direction</h2></summary>

In this section, let’s think only laterally (the refinement in the Z direction comes later). 
We are trying to find out how many grids do we want. From the previous part we have the minimum grid sizes. 

1. And for each casing we have its ID , therefore, the number of grids for each casing we should have round up casing_ID/min_grd_size.
    ```python
    no_grd_surf_case = math.ceil (surf_Casing_ID/min_grd_size)
    no_grd_cond_case = math.ceil (cond_Casing_ID/min_grd_size)
    no_grd_prod_case = math.ceil (prod_Casing_ID/min_grd_size)
    ```
    
2. Also, we do have the thickness of the cement bond from the previous part, accordingly we can assume how  many grids do we want that for this. 
    ```python
    no_grd_surf_bond = math.ceil (surf_bond/min_grd_size)
    no_grd_cond_bond = math.ceil (cond_bond/min_grd_size)
    no_grd_prod_bond = math.ceil  (prod_bond/min_grd_size)
    ```
    
3. Then, per casing the total number of girds = number of girds to mimic the pipe + 2 x number of grids for the cement bond. Thinking laterally we need one to the right and one to the left. 
    ```python
    no_latral_grd_surf = no_grd_surf_case + no_grd_surf_bond*2
    no_latral_grd_cond = no_grd_cond_case + no_grd_cond_bond*2
    no_latral_grd_prod = no_grd_prod_case + no_grd_prod_bond*2
    ```

4. We want our LGR to fine grids laterally to cover all of the casings from widest to narrowest one. So, we will have to find the max value of each calculated number of grids per casing.
    ```python
    no_latral_fine_grd = max (no_latral_grd_prod,no_latral_grd_cond,no_latral_grd_surf)
    ```



5. Now we create a list with only the smallest grids.   
    ```python
    LGR_size_fine_grd = [min_grd_size]*no_latral_fine_grd
    ```


6. we add 10 times and then 100 times larges grids around the fine ones, now we are moving from inner most grids to the outer most grids.

    ```python
    LGR_size_no_outer =[min_grd_size*100] +[min_grd_size*10] +[min_grd_size]*no_latral_fine_grd+ [min_grd_size*10] +[min_grd_size*100]
    ```

7. What is left behind divided by 2 could be placed in the outermost grids. 
    ```python
    LGR_sizes_xy = [(main_grd_dx- sum (LGR_size_no_outer))/2]+ LGR_size_no_outer + [(main_grd_dx- sum (LGR_size_no_outer))/2]
    ```
All of the well elements are located in the 0.05 m area, all the rest of grids is to fill-up the rest of the area gradually. 

![LGR result](https://github.com/equinor/GaP/blob/main/LGR_GaP.png)
   
</details>

<details><summary><h2>Number of fine grids in the whole LGR in Z direction</h2></summary>

After setting up the up the number of grids in X-Y direction, we should set up number of grids in Z direction.

The idea is not to refine too much in the Z direction, mainly because of the computational time. 

However, it we don’t refine at all, the resolution in the Z direction in overburden might be too low. Then, the it would limit our access to a specific depth, when we want to set start and end depth of each casing. 
First we define how many grids do we want? 
- We hard-coded the number of grids in overburden to 10 per layer. i.e. each layer is divided into the 10 layers. For example, if we have 10 grids in the overburden, in the LGR we will have 100 layers. 
- Then we convert the list to numpy array to be able to do calculation on that.
- Later we will see the LGR_num_z can directly be used in the NZFIN keyword under CARFIN 
   
   ```python
   LGR_numb_z = (no_of_layers_in_OB-main_grd_min_k + 1) * [10] + (main_grd_max_k -no_of_layers_in_OB) *[1] 
   LGR_numb_z = np.array (LGR_numb_z)
   ```



