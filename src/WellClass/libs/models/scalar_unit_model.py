
from pydantic import BaseModel

class ScalarUnitModel(BaseModel):
    """ model for scalar/unit pair
        Args:
            value (float): float value
            unit (str): name of metric unit
    """
    value: float|int
    unit: str
