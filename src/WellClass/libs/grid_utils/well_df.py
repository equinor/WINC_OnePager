
import pandas as pd

class WellDataFrame:

    def __init__(self, my_well):
        """ This is a container class
        """

        # # Dataframes for drilling, casings, borehole and geological tops

        self.drilling_df = pd.DataFrame(my_well.drilling)
        self.casings_df = pd.DataFrame(my_well.casings)
        self.borehole_df = pd.DataFrame(my_well.borehole)
        self.annulus_df = pd.DataFrame(my_well.annulus)
        self.geology_df = pd.DataFrame(my_well.geology)
        
        # and for the barriers
        self.barriers_df = pd.DataFrame(my_well.barriers)
        self.barriers_mod_df = pd.DataFrame(my_well.barriers_mod)
