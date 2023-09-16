
# handle type hints problem for python version < 3.10
from typing import Union, List

from pydantic import BaseModel

from .model_pipe_cement import PipeCementModel, ElemModel
from .model_defaults import DefaultModel
from .model_file import FileModel

class GaPSpec(BaseModel):
    """ specification for GaP
        Args:
            defaults (DefaultModel): some default values
            sim_case (FileModel): input file information for simulation case
            lgr_out (FileModel): output file information for LGR .grdecl
            casings (list[PipeCementModel]): information about casings
            barriers (list[ElemModel]): information about barriers, optional
    """
    defaults: DefaultModel

    sim_case: FileModel
    lgr_out: FileModel
    
    casings: List[PipeCementModel]
    barriers: Union[List[ElemModel], None] = None

class MetaDataModel(BaseModel):
    """ meta data
        Args:
            namespace (str): name space
            name (str): can use it for project name
            author (str): who made this yaml file
    """
    namespace: str = 'screen'
    name: Union[str, None] = None
    author: Union[str, None] = None

class GaPModel(BaseModel):
    """ contains all necessary parameters
        Args:
            apiVersion (str): current version of this yaml format
            kind (str): for GaP
            metadata (MetaDataModel): miscelaneous data
            spec (GapSpec): specification to GaP
    """
    apiVersion: str = 'gap/v0.1'
    kind: str = 'GaP'
    metadata: Union[MetaDataModel, None] = None
    spec: GaPSpec
