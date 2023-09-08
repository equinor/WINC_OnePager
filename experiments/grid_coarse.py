
import functools

import numpy as np

from ecl.eclfile import EclFile
from ecl.eclfile import EclInitFile, EclRestartFile
from ecl.grid import EclGrid
#import rips 

class GridCoarse:

    def __init__(self, simcase):
        """ initialize the grid
        """

        #Get grid dimensions and coordinates
        grid = EclGrid(simcase + ".EGRID") 
        #init = EclGrid(simcase + ".INIT") 
        init = EclInitFile(grid, simcase + ".INIT")
        # restart file
        # rst = EclRestartFile(grid, simcase + ".UNRST")

        # ## The grid dimensions
        self.NX, self.NY, self.NZ, _ = grid.get_dims()

        # # Store INIT parameters into a Pandas Dataframe: grid_init

        self.grid_init = grid.export_index()

        # Static properties Dataframe
        for key in init.keys():
                try:
                        self.grid_init[key] = init[key][0].numpy_view()
                except:
                        continue

    def set_cell_coords(self):
        """ Create cell coordinate X, Y, Z
        """

        grid_init = self.grid_init

        # generate cell coordinates by shifting half cell size
        xcoord = (grid_init.query("j==0&k==0").DX.cumsum() - grid_init.query("j==0&k==0").DX/2).values
        ycoord = (grid_init.query("i==0&k==0").DY.cumsum() - grid_init.query("i==0&k==0").DY/2).values
        zcoord = (grid_init.query("i==0&j==0").DZ.cumsum() - grid_init.query("i==0&j==0").DZ/2).values

        # TODO(hzh): a bug?
        # map_X = dict(zip(grid_init.query("j==0&j==0")['i'], xcoord))
        map_X = dict(zip(grid_init.query("j==0&k==0")['i'], xcoord))
        map_Y = dict(zip(grid_init.query("i==0&k==0")['j'], ycoord))
        map_Z = dict(zip(grid_init.query("i==0&j==0")['k'], zcoord))

        # save cell coordinates to DataFrame
        grid_init['X'] = grid_init['i'].map(map_X)
        grid_init['Y'] = grid_init['j'].map(map_Y)
        grid_init['Z'] = grid_init['k'].map(map_Z)

    def set_grid_info(self, Ali_way):
        """ Grid information for coarse grid
        """

        grid_init = self.grid_init
        
        # Retrieve coarse x-y grid indexes where LGR will be placed
        # i.e., center grid
        self.main_grd_i = grid_init.i.max()//2
        self.main_grd_j = grid_init.j.max()//2

        # Rettrieve min and max k-index for column where LGR will be placed
        self.main_grd_min_k = grid_init.k.min()
        self.main_grd_max_k = grid_init.k.max()

        # Retrieve coarse cell sizes
        self.main_grd_dx = grid_init.query('i == @main_grd_i & j == @main_grd_j & k == k.min()')['DX'].iloc[0]
        self.main_grd_dy = grid_init.query('i == @main_grd_i & j == @main_grd_j & k == k.min()')['DY'].iloc[0]


        # Retrieve all DZ in coarse grid, not used
        main_DZ = grid_init.query('i == @main_grd_i & j == @main_grd_j')['DZ'].values

        main_DEPTH = grid_init.query('i == @main_grd_i & j == @main_grd_j')['DEPTH'].values

        # depth where LGR starts
        self.ref_depth = 0
        if Ali_way:
            self.ref_depth = main_DEPTH[0] - 0.5*main_DZ[0]

        #Retrieve number of cells representing water column and overburden
        self.no_of_layers_in_OB    = grid_init.query('i==@main_grd_i & j == @main_grd_j & DZ >  10')['DZ'].shape[0]
        self.no_of_layers_below_OB = grid_init.query('i==@main_grd_i & j == @main_grd_j & DZ <= 10')['DZ'].shape[0]

        print('====>', no_of_layers_in_OB, no_of_layers_below_OB)        

    def set_DZ_rsrv_ovb(self):
        """ DZs for rsrv and ovb
        """
        grid_init = self.grid_init

        # the dz value to distinguish zones between reservoir and ovb
        dz0 = 10

        # center grid
        main_grd_i = self.main_grd_i
        main_grd_j = self.main_grd_j

        # 3.1 DZs for reservoir
        DZ_rsrv = grid_init.query('i==@main_grd_i & j == @main_grd_j & DZ <= @dz0')['DZ'].values

        # 3.2 DZs for coarse grid
        DZ_ovb_coarse = grid_init.query('i==@main_grd_i & j == @main_grd_j & DZ > @dz0')['DZ'].values

        return DZ_rsrv, DZ_ovb_coarse
    
    @functools.cached_property
    def main_grd_i(self):
        return self.grid_init.i.max()//2

    @functools.cached_property   
    def main_grd_j(self):
        return self.grid_init.j.max()//2
    
    def plot_raw(self, ):
        """ Plot well sketch and 2D slice of the permeability, at coarse grid
        """
        
        # for convenience
        grid_init = self.grid_init

        # middle indices of x and y
        mid_i = grid_init.i.max()//2
        mid_j = grid_init.j.max()//2

        # generate grid coordinates for plotting

        # grid coordinates
        xcorn  = (grid_init.query("j==0&k==0").DX.cumsum()).values
        ycorn  = (grid_init.query("i==0&k==0").DY.cumsum()).values
        zcorn  = (grid_init.query("i==0&j==0").DZ.cumsum()).values

        # add origin coordinates
        xcorn = np.append(0, xcorn)
        ycorn = np.append(0, ycorn)
        zcorn = np.append(0, zcorn)

        # shift grid coordinates half-length in x-y directions, i.e., [0, 3900] => [-1900, 2100]
        # but not in z direction
        xcorn -= xcoord[mid_i]
        ycorn -= ycoord[mid_j]

        # extract 2D xz slice at middle of y
        XZ_slice = grid_init.query('j==@mid_j')

        # extract permeability
        Z = XZ_slice.PERMX.values.reshape(NZ, NX)

        # plot x-z slice
        plot_well_perm(my_well, x=xcorn, y=zcorn, Z=Z, on_coarse=True)

