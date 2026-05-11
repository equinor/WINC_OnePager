"""Tests for utility modules: fraction_float, compute_intersection, csv_parser, yaml_parser."""

import json
import pathlib

import numpy as np
import pytest

from src.WellClass.libs.utils.compute_intersection import compute_intersection
from src.WellClass.libs.utils.csv_parser import csv_parser
from src.WellClass.libs.utils.fraction_float import float_to_fraction_inches, fraction_float
from src.WellClass.libs.utils.yaml_parser import yaml_parser

# ── fraction_float ──────────────────────────────────────────────────────


class TestFractionFloat:
    def test_integer_string(self):
        assert fraction_float("12") == 12

    def test_fraction_string(self):
        assert fraction_float("12 1/4") == 12.25

    def test_fraction_three_eighths(self):
        assert fraction_float("9 3/8") == 9.375

    def test_fraction_five_eighths(self):
        assert fraction_float("9 5/8") == 9.625

    def test_fraction_half(self):
        assert fraction_float("17 1/2") == 17.5


class TestFloatToFractionInches:
    def test_whole_number(self):
        assert float_to_fraction_inches(12.0) == '12"'

    def test_quarter(self):
        result = float_to_fraction_inches(12.25)
        assert result == '12 1/4"'

    def test_half(self):
        result = float_to_fraction_inches(17.5)
        assert result == '17 1/2"'

    def test_roundtrip(self):
        """fraction_float and float_to_fraction_inches should be inverse operations."""
        original = "9 5/8"
        as_float = fraction_float(original)
        back = float_to_fraction_inches(as_float)
        # strip the trailing "
        assert back == f'{original}"'


# ── compute_intersection ────────────────────────────────────────────────


class TestComputeIntersection:
    def test_simple_crossing(self):
        x = np.array([0.0, 1.0, 2.0, 3.0])
        y1 = np.array([0.0, 1.0, 2.0, 3.0])
        y2 = np.array([3.0, 2.0, 1.0, 0.0])
        xi, yi = compute_intersection(x, y1, y2)
        assert pytest.approx(xi, abs=0.01) == 1.5
        assert pytest.approx(yi, abs=0.01) == 1.5

    def test_no_intersection(self):
        x = np.array([0.0, 1.0, 2.0])
        y1 = np.array([1.0, 2.0, 3.0])
        y2 = np.array([5.0, 6.0, 7.0])
        xi, yi = compute_intersection(x, y1, y2)
        assert np.isnan(xi)
        assert np.isnan(yi)

    def test_multiple_intersections_returns_largest_x(self):
        x = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
        y1 = np.array([0.0, 2.0, 0.0, 2.0, 0.0])
        y2 = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
        xi, yi = compute_intersection(x, y1, y2)
        # Should return the rightmost intersection
        assert xi > 2.0

    def test_handles_nan_in_input(self):
        x = np.array([0.0, 1.0, np.nan, 3.0, 4.0])
        y1 = np.array([0.0, 1.0, np.nan, 3.0, 4.0])
        y2 = np.array([4.0, 3.0, np.nan, 1.0, 0.0])
        xi, yi = compute_intersection(x, y1, y2)
        assert np.isfinite(xi)
        assert np.isfinite(yi)


# ── csv_parser ──────────────────────────────────────────────────────────


class TestCsvParser:
    @pytest.fixture
    def csv_path(self):
        return pathlib.Path("test_data/examples/wildcat/GaP_input_Wildcat_v3.csv")

    def test_returns_dict(self, csv_path):
        result = csv_parser(csv_path)
        assert isinstance(result, dict)

    def test_expected_keys(self, csv_path):
        result = csv_parser(csv_path)
        expected_keys = {"well_header", "drilling", "casing_cement", "barriers", "geology"}
        assert expected_keys.issubset(result.keys())

    def test_well_header_has_well_name(self, csv_path):
        result = csv_parser(csv_path)
        assert "well_name" in result["well_header"]

    def test_drilling_has_diameter(self, csv_path):
        result = csv_parser(csv_path)
        assert "diameter_in" in result["drilling"]


# ── yaml_parser ─────────────────────────────────────────────────────────


class TestYamlParser:
    @pytest.fixture
    def yaml_path(self):
        return pathlib.Path("test_data/examples/simple_well/Simple_well.yaml")

    @pytest.fixture
    def wildcat_yaml_path(self):
        return pathlib.Path("test_data/examples/wildcat/wildcat.yaml")

    def test_returns_well_model(self, yaml_path):
        result = yaml_parser(yaml_path)
        assert hasattr(result, "spec")

    def test_well_name_parsed(self, yaml_path):
        result = yaml_parser(yaml_path)
        spec = json.loads(result.spec.model_dump_json())
        assert spec["well_header"]["well_name"] == "WellA"

    def test_drilling_sections_parsed(self, yaml_path):
        result = yaml_parser(yaml_path)
        spec = json.loads(result.spec.model_dump_json())
        assert len(spec["drilling"]) == 3

    def test_barriers_parsed(self, yaml_path):
        result = yaml_parser(yaml_path)
        spec = json.loads(result.spec.model_dump_json())
        assert len(spec["barriers"]) == 2

    def test_wildcat_yaml(self, wildcat_yaml_path):
        result = yaml_parser(wildcat_yaml_path)
        spec = json.loads(result.spec.model_dump_json())
        assert "well_header" in spec
        assert spec["co2_datum"] == 2370
