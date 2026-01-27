from src.WellClass.libs.well_class import Well
from tests.conftest import example_well_dict

def test_well_header(well_class_fixture,
                     well_class_dict_fixture) -> None:
    """ test well_header
        Args:
            well_class_fixture: fixture for WellClass
            well_class_dict_fixture: fixture for WellClass example
    """

    info    = well_class_fixture.header
    info_gt = well_class_dict_fixture['header']

    assert info['well_name'] == info_gt['well_name']
    assert info['well_td_rkb'] == info_gt['well_td_rkb']

    # compare the whole dictionary
    assert info == info_gt

def test_drilling(well_class_fixture,
                  well_class_dict_fixture) -> None:
    """ test drilling
        Args:
            well_class_fixture: fixture for WellClass
            well_class_dict_fixture: fixture for WellClass example
    """

    info = well_class_fixture.drilling
    info_gt = well_class_dict_fixture['drilling']

    # section index
    sec_index = 1

    assert info['diameter_in'][sec_index] == info_gt['diameter_in'][sec_index]
    assert info['top_rkb'][sec_index] == info_gt['top_rkb'][sec_index]
    assert info['bottom_rkb'][sec_index] == info_gt['bottom_rkb'][sec_index]

    # compare the whole dictionary
    assert info == info_gt

def test_casings(well_class_fixture,
                 well_class_dict_fixture) -> None:
    """ test casings
        Args:
            well_class_fixture: fixture for WellClass
            well_class_dict_fixture: fixture for WellClass example
    """

    info = well_class_fixture.casings
    info_gt = well_class_dict_fixture['casings']

    # section index
    sec_index = 1

    assert info['diameter_in'][sec_index] == info_gt['diameter_in'][sec_index]
    assert info['top_rkb'][sec_index] == info_gt['top_rkb'][sec_index]
    assert info['bottom_rkb'][sec_index] == info_gt['bottom_rkb'][sec_index]

    # compare the whole dictionary
    assert info == info_gt

def test_barriers(well_class_fixture,
                  well_class_dict_fixture) -> None:
    """ test barriers
        Args:
            well_class_fixture: fixture for WellClass
            well_class_dict_fixture: fixture for WellClass example
    """

    info = well_class_fixture.barriers
    info_gt = well_class_dict_fixture['barriers']

    # section index
    sec_index = 1

    assert info['barrier_name'][sec_index] == info_gt['barrier_name'][sec_index]
    assert info['top_rkb'][sec_index] == info_gt['top_rkb'][sec_index]
    assert info['bottom_rkb'][sec_index] == info_gt['bottom_rkb'][sec_index]

    # compare the whole dictionary
    assert info == info_gt
    

def test_well_instantiation_from_dict(well_class_fixture,
                                      well_class_dict_fixture):

    assert well_class_fixture.header['well_name'] == well_class_dict_fixture['header']['well_name']
    assert isinstance(well_class_fixture.casings, dict)
    assert well_class_fixture.casings['diameter_in']  # Check key exists


def test_well_instantiation_from_yaml(well_class_fixture):
    assert well_class_fixture.header['well_name']
    assert well_class_fixture.barriers  # Should not be empty