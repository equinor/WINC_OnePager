'''
How to initialize:

```
import well_class

INDATA          = <path to a csv-file with the well-data>

#Reads the csv-file and organize the data into a dict of dataframes

well_df         = well_class.csv_parser(INDATA)                          
```

Then the class is initialized with a lot of explicit calls. (Bad structure - should been done in one go: ```mywell = Well(INDATA)```)
```
mywell          = well_class.Well(
                       header       = well_df['well_header'],
                       reservoir_P  = well_df['reservoir_pressure'],
                       drilling     = well_df['drilling'],
                       casings      = well_df['casing_cement'],
                       barriers     = well_df['barriers'],
                       geology      = well_df['geology'],
                       main_barrier = well_df['main_barrier'],
                       barrier_perm = well_df['barrier_permeability'],
                       co2_datum    = well_df['co2_datum']
                   )
```

Now additional functionalities that can be explicitely called are 
```
   .plot_pt()

   .plot_pressure()  + plt.show()

   .plot_sketch()    + plt.show()
```

'''

from dataclasses import dataclass
import json

from ..well_computed import (
    compute_borehole,
    compute_cement_bond,
    compute_annulus,
    compute_barriers_diam,
    get_barriers_names,
    get_barrier_height_and_depth,
    get_barrier_radius,
)

from .well_raw import WellRaw

@dataclass              # @dataclass(kw_only=True)
class Well(WellRaw):
    """ This contains not only the basic well information but also its computed information.

        Args:
            borehold (dict): for borehole information
            cement_bond (dict): contains information about cement bond
            annulus (dict): gap between casing and openhole
            barriers_mod (dict): extra information about barriers
            barriers_names (dict): reorgainze barrier names
    """
    borehole      : dict = None
    cement_bond   : dict = None
    annulus       : dict = None
    barriers_mod  : dict = None
    barriers_names: dict = None

    def __post_init__(self):

        super().__post_init__()

        self._compute_well()

    def _compute_well(self):
        """ compute extra well information
        """
        
        self.borehole = compute_borehole(self.casings, self.drilling)
        self.cement_bond = compute_cement_bond(self.casings, self.drilling)
        self.annulus= compute_annulus(self.casings, self.drilling)
        
        self.barriers_mod = compute_barriers_diam(self.barriers, self.borehole)
        self.barriers_names = get_barriers_names(self.barriers_mod)

    def compute_barrier_props(self, barrier_name: str) -> dict:
        """ Compute barrier geometrical information

            Args:
               barrier_name (str): barrier name 
        """

        # for convenience
        barriers_mod = self.barriers_mod
        barriers_names = self.barriers_names

        # properties
        barrier_props = {}

        # height/depth
        barrier_h_d = get_barrier_height_and_depth(barriers_mod, barriers_names, barrier_name)
        barrier_props.update(barrier_h_d)
        
        # radius
        barrier_r = get_barrier_radius(barriers_mod, barriers_names, barrier_name)
        barrier_props.update(barrier_r)

        return barrier_props

    @property
    def to_json(self):
        return json.dumps(self.__dict__, indent=4)
