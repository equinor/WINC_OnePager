
# handle type hints problem for python version < 3.10
from typing import Union

from pydantic import BaseModel

class ScalarUnitModel(BaseModel):
    """ model for scalar/unit pair
        Args:
            value (float): float value
            unit (str): name of metric unit
    """
    value: Union[float, int]
    unit: str
