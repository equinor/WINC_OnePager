
import pytest

from src.WellClass.libs.grid_utils import (
    WellDataFrame,
)

@pytest.fixture(scope='session')
def well_dataframe_fixture(well_class_fixture):
    """ Well DataFrame fixture
        Args:
            well_class_fixture: fixture for WellClass
    """
    well_df = WellDataFrame(well_class_fixture)

    return well_df
