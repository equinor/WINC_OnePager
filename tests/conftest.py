""" define some fixtures
"""
import pathlib
import json
import pickle

import pytest

from src.WellClass.libs.utils import (
    csv_parser,
    yaml_parser,
)

from src.WellClass.libs.well_class import Well

# a test example: cosmo
example = {
    'sim_path': './test_data/examples/cosmo',
    'well_config': 'cosmo.yaml',
    'sim_case': 'TEMP-0.EGRID'
}

@pytest.fixture(scope='session')
def well_class_fixture():
    """ fixture for WellClass
    """

    # the paths
    sim_path = pathlib.Path(example['sim_path'])
    well_name = sim_path/example['well_config']

    # extract suffix
    suffix = well_name.suffix
    # .yaml or .csv?
    use_yaml = suffix in ['.yaml', '.yml']

    if use_yaml:
        
        # # pydantic model
        well_model = yaml_parser(well_name)
        well_csv = json.loads(well_model.spec.model_dump_json())
    else:

        # load the well information
        well_csv = csv_parser(well_name)

    # instantiate class
    my_well = Well(header       = well_csv['well_header'], 
                   drilling     = well_csv['drilling'],
                   casings      = well_csv['casing_cement'],
                   geology      = well_csv['geology'],
                   barriers     = well_csv['barriers'], 
                   barrier_perm = well_csv['barrier_permeability'],
                   co2_datum    = well_csv['co2_datum'],
                )
    
    return my_well

@pytest.fixture(scope='session')
def well_class_dict_fixture():
    """ fixture for loading WellClass .pkl file for unit testing
    """
    # PYTEST folder
    PYTEST_FOLDER = 'pytest'

    # the paths
    pytest_path = pathlib.Path(example['sim_path'])/PYTEST_FOLDER
    well_name = pytest_path/example['well_config']

    # extract suffix
    suffix = well_name.suffix
    # .yaml or .csv?
    use_yaml = suffix in ['.yaml', '.yml']

    # only handle .yaml file
    assert use_yaml is True

    # convert .yaml file suffix to .pkl
    well_name_tuple = str(well_name).split('.')
    well_name_pkl = well_name_tuple[0]+'_well_pytest.pkl'

    # load .pkl file
    with open(well_name_pkl, 'rb') as f:
        my_well = f.read()

    # convert it to dictionary
    my_well_dict = pickle.loads(my_well)

    return my_well_dict

