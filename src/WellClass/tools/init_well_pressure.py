# WellClass
# def init_well_pressure():
#     # input configuration file name,
#     # for example, './test_data/examples/smeaheia_v1.yaml'
#     well_name = Path(args.config_file)
#     # input mixture index.
#     mixture_index = args.mixture_index
#     print(f'{mixture_index=}')
#     mixtures_list = ['pure_co2', 'mixture1', 'mixture2']
#     mixture = mixtures_list[mixture_index]
#     pvt_files = files('src.WellClass.libs.pvt.pvt_constants')
#     pvt_path = pvt_files / mixture
#     #depth values to evaluyate pressure at Shmin
#     maxP_values = args.max_pressure_depths
#     # extract suffix
#     suffix = well_name.suffix
#     # .yaml or .csv?
#     use_yaml = suffix in ['.yaml', '.yml']
#     ############ 2. Load well configuration file ###############
#     if use_yaml:
#         # # pydantic model
#         well_model = yaml_parser(well_name)
#         well_csv = json.loads(well_model.spec.model_dump_json())
#     else:
#         # load the well information
#         well_csv = csv_parser(well_name)
#     ########### 3. build Well and pressure classes ######################
#     # 3.1 build well class
#     my_well = Well( header       = well_csv['well_header'],
#                     co2_datum    = well_csv['co2_datum'],
#             )
#     # 3.2 build pressure class
#     my_pressure = Pressure( header      = well_csv['well_header'],
#                             reservoir_P = well_csv['reservoir_pressure'],
#                             co2_datum   = well_csv['co2_datum'],
#                             max_pressure_pos = maxP_values,
#                             pvt_path    = pvt_path,)
import json
from pathlib import Path
from typing import Any

# parsing libraries
from ..libs.utils import (
    csv_parser,
    yaml_parser,
)
from ..libs.well_class import Well

# Pressure Class


def initialize_well(
    indata: str | None = None,
    header: dict[str, Any] | None = None,
    reservoir_P: dict[str, Any] | None = None,
    drilling: dict[str, Any] | None = None,
    casings: dict[str, Any] | None = None,
    barriers: dict[str, Any] | None = None,
    geology: dict[str, Any] | None = None,
    main_barrier: dict[str, Any] | None = None,
    barrier_perm: dict[str, Any] | None = None,
    co2_datum: float | None = None,
) -> Well:
    """
    Function to initialize the Well class. It can either accept dictionaries directly,
    or a YAML file path to load the dictionaries.

    :param indata: Path to a YAML file containing the well data
    :param header: Dictionary with well header data
    :param reservoir_P: Dictionary with reservoir pressure data
    :param drilling: Dictionary with drilling data
    :param casings: Dictionary with casing cement data
    :param barriers: Dictionary with barrier data
    :param geology: Dictionary with geology data
    :param main_barrier: Dictionary with main barrier data
    :param barrier_perm: Dictionary with barrier permeability data
    :param co2_datum: Co2 datum value
    :return: A Well object initialized with the provided data
    """
    if indata:
        well_path = Path(indata)
        # extract suffix
        suffix = well_path.suffix
        # .yaml or .csv?
        use_yaml = suffix in [".yaml", ".yml"]

        if use_yaml:
            # # pydantic model
            well_model = yaml_parser(well_path)
            well_csv = json.loads(well_model.spec.model_dump_json())
        else:
            # load the well information
            well_csv = csv_parser(well_path)

        header = well_csv["well_header"]
        drilling = well_csv["drilling"]
        casings = well_csv["casing_cement"]
        barriers = well_csv["barriers"]
        barrier_perm = well_csv["barriers_permeability"]
        geology = well_csv["geology"]
        co2_datum = well_csv["co2_datum"]

    # At this point, regardless of the input method, we should have all required dictionaries or None

    # Build Well Class
    mywell = Well(
        header=header, drilling=drilling, casings=casings, barriers=barriers, barrier_perm=barrier_perm, geology=geology, co2_datum=co2_datum
    )

    # If you had a separate class that needed the reservoir_P, main_barrier, etc., you would Initialize it here

    return mywell
