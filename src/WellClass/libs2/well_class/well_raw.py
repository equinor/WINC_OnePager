from dataclasses import dataclass
from typing import Any


@dataclass
class WellRaw:
    """
    Basic user input well information
    Args:
    header (dict): well header
    hole_casings (list of dict): list of hole and casing information
    survey (dict of list): survey information
    plugs (list of dict): list of plug information
    stratigraphy (list of dict): list of stratigraphy information
    """

    header: dict[str, Any] | None = None
    hole_casings: list[dict[str, Any]] | None = None
    survey: dict[str, list[float]] | None = None
    plugs: list[dict[str, Any]] | None = None
    stratigraphy: list[dict[str, Any]] | None = None
