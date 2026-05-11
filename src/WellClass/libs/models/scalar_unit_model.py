# handle type hints problem for python version < 3.10

from pydantic import BaseModel


class ScalarUnitModel(BaseModel):
    """
    model for scalar/unit pair
    Args:
        value (float): float value
        unit (str): name of metric unit
    """

    value: float | int
    unit: str
