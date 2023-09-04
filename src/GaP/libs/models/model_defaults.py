
from pydantic import BaseModel

class DefaultModel(BaseModel):
    """ contains some defaults
        Args:
            mindz_ob: minimum dz value to be considered to be overburden
    """
    mindz_ob: float = 10.0
