from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import wellpathpy as wp
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from scipy import constants as const

from src.WellClass.libs2.plotting.plot_sketch import plot_sketch
from src.WellClass.libs2.well_class import Well
from src.WellClass.libs2.well_class.well_validation import (
    split_hole_casings,
    verify_hole_casings,
    verify_plugs,
    verify_stratigraphy,
)
from src.WellClass.libs2.well_computed.annulus import compute_annulus
from src.WellClass.libs2.well_computed.borehole import compute_borehole
from src.WellClass.libs2.well_computed.plug_diameter import compute_plug_diameter
from src.WellClass.libs2.well_computed.plug_properties import compute_plug_properties
from src.WellClass.libs2.well_computed.well_path import (
    build_wellpath_object,
    md2tvd_interpolator,
)


@dataclass
class WellProcessed(Well):
    wellpath: wp.position_log | None = None
    md2tvd: Callable[[float], float] | None = None
    borehole: list[dict[str, Any]] | None = None
    cement_bond: list[dict[str, Any]] | None = None
    annulus: list[dict[str, Any]] | None = None
    processed_plugs: list[dict[str, Any]] | None = None
    plug_names: list[str] | None = None

    def __post_init__(self) -> None:
        # invoke parent post_init to compute wellpath and md2tvd
        super().__post_init__()

        self._check_header_units()

        self.wellpath = self._build_wellpath()
        self.md2tvd = self._md2tvd_interpolator()

        if self.inventory["hole_casings"]:
            self._process_hole_casings()

        if self.inventory["plugs"]:
            self._process_plugs()

        if self.inventory["stratigraphy"]:
            self._process_stratigraphy()

        self._compute_well()

    def _check_header_units(self) -> None:
        if not self.header:
            raise ValueError("Header is required to check units")

        if self.header.get("depth_reference_rkb_unit") == "ft":
            self.header["depth_reference_rkb_unit"] = "m"
            self.header["total_depth_rkb"] *= const.foot

        if self.header.get("ground_elevation_unit") == "ft":
            self.header["ground_elevation_unit"] = "m"
            self.header["ground_elevation"] *= const.foot

        if self.header.get("total_depth_rkb_unit") == "ft":
            self.header["total_depth_rkb_unit"] = "m"
            self.header["total_depth_rkb"] *= const.foot

    def _build_wellpath(self) -> wp.position_log:
        return build_wellpath_object(survey=self.survey, total_depth=self.header["total_depth_rkb"], survey_bool=self.inventory["survey"])

    def _md2tvd_interpolator(self) -> Callable[[float], float]:
        """Create an interpolator function to convert md to tvd"""
        return md2tvd_interpolator(self.wellpath, self.header["depth_reference_rkb"])

    def _process_intervals(self, intervals: list[dict], include_diameter: bool = False) -> list[dict]:
        processed = []
        for interval in intervals:
            top_tvd_msl = self.md2tvd(interval["top_rkb"])
            bottom_tvd_msl = self.md2tvd(interval["bottom_rkb"])
            new_data = {
                **interval,
                "tvd_msl_top": top_tvd_msl,
                "tvd_msl_bottom": bottom_tvd_msl,
            }

            if include_diameter:
                new_data["diameter_m"] = interval["diameter_in"] * const.inch

            processed.append(new_data)

        return processed

    def _process_hole_casings(self) -> None:
        if not self.hole_casings:
            return
        processed = self._process_intervals(self.hole_casings, include_diameter=True)
        verify_hole_casings(processed)
        self.hole_casings = processed

    def _process_plugs(self) -> None:
        if not self.plugs:
            return
        processed = self._process_intervals(self.plugs)
        verify_plugs(processed, self.header["ground_elevation"])
        self.plugs = processed

    def _process_stratigraphy(self) -> None:
        if not self.stratigraphy:
            return

        processed = self._process_intervals(self.stratigraphy)
        verify_stratigraphy(processed, self.header["ground_elevation"])
        self.stratigraphy = processed

    def _compute_well(self) -> None:
        splitted_hole_casing = split_hole_casings(self.hole_casings)

        casings = splitted_hole_casing["casing"]
        holes = splitted_hole_casing["holes"]
        cement_bond = splitted_hole_casing["casing_cement"]

        self.borehole = compute_borehole(holes, casings, self.md2tvd)
        self.cement_bond = compute_annulus(holes=holes, casings=casings, cement_bond=cement_bond, md2tvd=self.md2tvd, solve_cement_bond=True)
        self.annulus = compute_annulus(holes=holes, casings=casings, cement_bond=cement_bond, md2tvd=self.md2tvd)

        if self.inventory["plugs"]:
            self.processed_plugs = compute_plug_diameter(self.plugs, self.borehole)
            self.plug_names = [plug["name"] for plug in self.plugs]

    def _compute_plug_properties(self, eval_plug: str) -> dict[str, Any]:
        return compute_plug_properties(self.processed_plugs, self.plug_names, eval_plug)

    def plot_sketch(self, ax: Axes | None = None, **kwargs: str | float | bool) -> tuple[Figure, Axes]:
        fig, ax = plot_sketch(self, ax=ax, **kwargs)

        return fig, ax
