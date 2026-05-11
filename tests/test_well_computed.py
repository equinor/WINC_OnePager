"""Tests for well_computed modules: borehole, annulus, cement_bond, barriers_mod, barriers_names, barrier_props."""

import json
import pathlib

import pandas as pd
import pytest

from src.WellClass.libs.utils.yaml_parser import yaml_parser
from src.WellClass.libs.well_class.well_class import Well
from src.WellClass.libs.well_computed.barrier_props import (
    compute_barrier_props,
)


@pytest.fixture(scope="module")
def simple_well():
    yaml_path = pathlib.Path("test_data/examples/wildcat/wildcat.yaml")
    well_model = yaml_parser(yaml_path)
    well_csv = json.loads(well_model.spec.model_dump_json())
    return Well(
        header=well_csv["well_header"],
        drilling=well_csv["drilling"],
        casings=well_csv["casing_cement"],
        geology=well_csv["geology"],
        barriers=well_csv["barriers"],
        barrier_perm=well_csv.get("barrier_permeability"),
        co2_datum=well_csv["co2_datum"],
    )


# ── compute_borehole ───────────────────────────────────────────────────


class TestComputeBorehole:
    def test_returns_dict(self, simple_well):
        borehole = simple_well.borehole
        assert isinstance(borehole, dict)

    def test_has_expected_keys(self, simple_well):
        borehole = simple_well.borehole
        assert "top_msl" in borehole
        assert "bottom_msl" in borehole
        assert "diameter_m" in borehole

    def test_diameters_are_decreasing(self, simple_well):
        borehole = simple_well.borehole
        diams = list(borehole["diameter_m"].values())
        # Each subsequent section should be smaller or equal
        for i in range(1, len(diams)):
            assert diams[i] <= diams[i - 1]

    def test_depths_are_continuous(self, simple_well):
        borehole = simple_well.borehole
        borehole_df = pd.DataFrame(borehole)
        # Bottom of each section should match top of next (except last)
        for i in range(len(borehole_df) - 1):
            assert borehole_df.loc[i, "bottom_msl"] == pytest.approx(borehole_df.loc[i + 1, "top_msl"])


# ── compute_annulus ─────────────────────────────────────────────────────


class TestComputeAnnulus:
    def test_returns_dict(self, simple_well):
        annulus = simple_well.annulus
        assert isinstance(annulus, dict)

    def test_has_expected_keys(self, simple_well):
        annulus = simple_well.annulus
        for key in ["ann_od_m", "ann_id_m", "top_msl", "bottom_msl", "A_i", "A_o", "thick_m"]:
            assert key in annulus

    def test_outer_larger_than_inner(self, simple_well):
        annulus = simple_well.annulus
        for idx in annulus["ann_od_m"]:
            assert annulus["ann_od_m"][idx] > annulus["ann_id_m"][idx]

    def test_thickness_is_positive(self, simple_well):
        annulus = simple_well.annulus
        for idx in annulus["thick_m"]:
            assert annulus["thick_m"][idx] > 0


# ── compute_cement_bond ─────────────────────────────────────────────────


class TestComputeCementBond:
    def test_returns_dict(self, simple_well):
        cb = simple_well.cement_bond
        assert isinstance(cb, dict)

    def test_has_expected_keys(self, simple_well):
        cb = simple_well.cement_bond
        for key in ["top_msl", "bottom_msl", "id_m", "od_m"]:
            assert key in cb

    def test_top_above_bottom(self, simple_well):
        cb = simple_well.cement_bond
        for idx in cb["top_msl"]:
            assert cb["top_msl"][idx] <= cb["bottom_msl"][idx]


# ── compute_barriers_diam ───────────────────────────────────────────────


class TestComputeBarriersDiam:
    def test_returns_dict(self, simple_well):
        barriers_mod = simple_well.barriers_mod
        assert isinstance(barriers_mod, dict)

    def test_has_expected_keys(self, simple_well):
        barriers_mod = simple_well.barriers_mod
        for key in ["b_name", "top_msl", "bottom_msl", "diameter_m"]:
            assert key in barriers_mod

    def test_barrier_names_present(self, simple_well):
        barriers_mod = simple_well.barriers_mod
        names = list(barriers_mod["b_name"].values())
        assert len(names) > 0


# ── get_barriers_names ──────────────────────────────────────────────────


class TestGetBarriersNames:
    def test_returns_dict(self, simple_well):
        barriers_names = simple_well.barriers_names
        assert isinstance(barriers_names, dict)

    def test_groups_sections(self, simple_well):
        barriers_names = simple_well.barriers_names
        # Each key should map to a list of section keys
        for name, sections in barriers_names.items():
            assert isinstance(sections, list)
            assert len(sections) >= 1


# ── barrier_props ───────────────────────────────────────────────────────


class TestBarrierProps:
    def test_compute_barrier_props(self, simple_well):
        barriers_mod = simple_well.barriers_mod
        barriers_names = simple_well.barriers_names
        # Pick the first barrier
        barrier_name = next(iter(barriers_names.keys()))
        props = compute_barrier_props(barriers_mod, barriers_names, barrier_name)
        assert "height" in props
        assert "top" in props
        assert "bottom" in props
        assert "radius" in props

    def test_height_is_positive(self, simple_well):
        barriers_mod = simple_well.barriers_mod
        barriers_names = simple_well.barriers_names
        for barrier_name in barriers_names:
            props = compute_barrier_props(barriers_mod, barriers_names, barrier_name)
            assert props["height"] > 0

    def test_radius_is_positive(self, simple_well):
        barriers_mod = simple_well.barriers_mod
        barriers_names = simple_well.barriers_names
        for barrier_name in barriers_names:
            props = compute_barrier_props(barriers_mod, barriers_names, barrier_name)
            assert props["radius"] > 0

    def test_top_above_bottom(self, simple_well):
        barriers_mod = simple_well.barriers_mod
        barriers_names = simple_well.barriers_names
        for barrier_name in barriers_names:
            props = compute_barrier_props(barriers_mod, barriers_names, barrier_name)
            assert props["top"] < props["bottom"]

    def test_height_matches_depth_span(self, simple_well):
        barriers_mod = simple_well.barriers_mod
        barriers_names = simple_well.barriers_names
        for barrier_name in barriers_names:
            props = compute_barrier_props(barriers_mod, barriers_names, barrier_name)
            assert props["height"] == pytest.approx(props["bottom"] - props["top"])


# ── Well class computed properties ──────────────────────────────────────


class TestWellComputedIntegration:
    def test_well_has_borehole(self, simple_well):
        assert simple_well.borehole is not None

    def test_well_has_annulus(self, simple_well):
        assert simple_well.annulus is not None

    def test_well_has_cement_bond(self, simple_well):
        assert simple_well.cement_bond is not None

    def test_well_has_barriers_mod(self, simple_well):
        assert simple_well.barriers_mod is not None

    def test_well_has_barriers_names(self, simple_well):
        assert simple_well.barriers_names is not None

    def test_well_inventory_all_true(self, simple_well):
        """Simple well has all sections, so all inventory flags should be True."""
        for key in ["drilling", "casings", "barriers", "geology"]:
            assert simple_well.inventory[key] is True

    def test_well_compute_barrier_props_method(self, simple_well):
        barrier_name = next(iter(simple_well.barriers_names.keys()))
        props = simple_well.compute_barrier_props(barrier_name)
        assert "height" in props
        assert "radius" in props

    def test_well_to_json(self, simple_well):
        json_str = simple_well.to_json
        assert isinstance(json_str, str)
        data = json.loads(json_str)
        assert "header" in data
