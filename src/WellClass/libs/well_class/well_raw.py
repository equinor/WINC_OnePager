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

# handle type hints problem for python version < 3.10
from typing import Union

from dataclasses import dataclass
import json

import pandas as pd
import scipy
import scipy.constants

from .well_raw_validation import (
    valid_drilling,
    valid_casings    
)

@dataclass              # @dataclass(kw_only=True)
class WellRaw:
    """ Basic user input well information

        Args:
            header (dict): well header
            drilling (dict): drilling well
            casings (dict): well casing
            barriers (dict): barrier information
            barrier_perm (dict): barrier permeabilities
            geology (dict): gelogical formation
            co2_datum (float): co2 datum value
    """
    header        : dict = None
    drilling      : dict = None
    casings       : dict = None
    barriers      : dict = None
    barrier_perm  : dict = None
    geology       : dict = None
    co2_datum     : Union[float, int] = None

    def __post_init__(self):
        """ compute basic well information
        """
        self._process_drilling()
        self._process_casings()
        self._process_barriers()
        self._process_geology()

    def _process_drilling(self):

        drilling_df = pd.DataFrame(self.drilling)

        drilling_df['diameter_m'] = drilling_df['diameter_in'] * scipy.constants.inch       #0.0254
        drilling_df['top_msl']    = drilling_df['top_rkb'] - self.header['well_rkb']
        drilling_df['bottom_msl'] = drilling_df['bottom_rkb'] - self.header['well_rkb']

        # validate drilling
        valid_drilling(drilling_df)

        self.drilling = drilling_df.to_dict()

    def _process_casings(self):

        casings_df = pd.DataFrame(self.casings)

        casings_df['diameter_m'] = casings_df['diameter_in'] * scipy.constants.inch         #0.0254
        casings_df['top_msl']    = casings_df['top_rkb']    - self.header['well_rkb']
        casings_df['bottom_msl'] = casings_df['bottom_rkb'] - self.header['well_rkb']
        casings_df['toc_msl']    = casings_df['toc_rkb']    - self.header['well_rkb']
        casings_df['boc_msl']    = casings_df['boc_rkb']    - self.header['well_rkb']

        # validate casings
        valid_casings(casings_df)

        self.casings = casings_df.to_dict()

    def _process_barriers(self):

        barriers_df = pd.DataFrame(self.barriers)

        barriers_df['top_msl']    = barriers_df['top_rkb']    - self.header['well_rkb']
        barriers_df['bottom_msl'] = barriers_df['bottom_rkb'] - self.header['well_rkb']

        # barriers_df.set_index('barrier_name', inplace=True)
        self.barriers = barriers_df.to_dict()
    
    def _process_geology(self):

        geology_df = pd.DataFrame(self.geology)

        geology_df = geology_df.dropna(how='all')
        geology_df = geology_df.reset_index(drop=True)

        geology_df['top_msl']  = geology_df['top_rkb'] - self.header['well_rkb']
        geology_df['base_msl'] = geology_df['top_msl'] - geology_df['top_msl'].diff(periods=-1)
        geology_df.loc[geology_df.index.max(), 'base_msl'] = self.header['well_td_rkb'] - self.header['well_rkb']

        self.geology = geology_df.to_dict()

    @property
    def to_json(self):
        return json.dumps(self.__dict__, indent=4)
