
from typing import Tuple
import numpy as np

from ecl.eclfile import EclFile, EclRestartFile
from ecl.eclfile import EclInitFile
from ecl.grid import EclGrid

from .extract_grid_utils import (
    extract_xz_corn_coords,
    extract_xz_prop_slice,
)

class GridLGR:

    def __init__(self, sim_case: str) -> None:

        #Get grid dimensions and coordinates
        grid = EclGrid(sim_case + ".EGRID") 
        #init = EclGrid(simcase + ".INIT") 
        init = EclInitFile(grid, sim_case + ".INIT")

        # the coarse grid dimensions
        self.NX, self.NY, self.NZ, _ = grid.get_dims()

        #Process init file
        lgr_name = grid.get_lgr(0).get_name()
        lgr_grid = grid.get_lgr(lgr_name)

        # # Store INIT parameters into a Pandas Dataframe

        lgr_index = lgr_grid.export_index()

        # Static properties Dataframe
        for key in init.keys():
            try:
                lgr_index[key] = init[key][1].numpy_view()
            except:
                continue

        # 
        self.lgr_index = lgr_index

        # compute DX and DY
        mid_i = lgr_index.i.max()//2
        mid_j = lgr_index.j.max()//2
 
        self.main_grd_dx = lgr_index.query("j==@mid_j&k==0").DX.sum()
        self.main_grd_dy = lgr_index.query("i==@mid_i&k==0").DY.sum()

    def extract_xz_corn_coords(self) -> Tuple[np.ndarray, np.ndarray]:
        """ generate xcorn and zcorn coordinates
        """

        # for convenience
        lgr_index = self.lgr_index

        # for shifting
        sDX = self.main_grd_dx/2
        sDY = self.main_grd_dy/2

        # generate grid coordinates for plotting
        xcorn, zcorn = extract_xz_corn_coords(lgr_index, sDX, sDY)

        return xcorn, zcorn 
    
    def extract_xz_slice(self) -> np.ndarray:
        """ generate x-z PERM slice
        """
        # for convenience
        lgr_index = self.lgr_index

        # extract permeability
        Z = extract_xz_prop_slice(lgr_index)

        return Z 
