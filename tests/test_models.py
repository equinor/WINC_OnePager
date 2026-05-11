"""Tests for Pydantic models: well_model, well_model_utils, scalar_unit_model."""

import pytest
from pydantic import ValidationError

from src.WellClass.libs.models.scalar_unit_model import ScalarUnitModel
from src.WellClass.libs.models.well_model import MetaDataModel, WellModel, WellPressureSpec, WellSpec
from src.WellClass.libs.models.well_model_utils import (
    BarrierModel,
    CasingCementModel,
    DrillingModel,
    GeologyModel,
    ReservoirPressureModel,
    WellHeaderModel,
)

# ── ScalarUnitModel ─────────────────────────────────────────────────────


class TestScalarUnitModel:
    def test_valid_creation(self):
        s = ScalarUnitModel(value=100.0, unit="m")
        assert s.value == 100.0
        assert s.unit == "m"

    def test_integer_value(self):
        s = ScalarUnitModel(value=42, unit="bar")
        assert s.value == 42


# ── WellHeaderModel ─────────────────────────────────────────────────────


class TestWellHeaderModel:
    def test_valid_header(self):
        h = WellHeaderModel(
            well_name="TestWell",
            well_rkb=25,
            sf_depth_msl=300,
            well_td_rkb=1500,
            sf_temp=4,
            geo_tgrad=40,
        )
        assert h.well_name == "TestWell"
        assert h.well_rkb == 25

    def test_missing_field_raises(self):
        with pytest.raises(ValidationError):
            WellHeaderModel(
                well_name="TestWell",
                well_rkb=25,
                # missing sf_depth_msl and others
            )


# ── DrillingModel ───────────────────────────────────────────────────────


class TestDrillingModel:
    def test_numeric_diameter(self):
        d = DrillingModel(top_rkb=100, bottom_rkb=500, diameter_in=12.25)
        assert d.diameter_in == 12.25

    def test_string_fraction_diameter(self):
        d = DrillingModel(top_rkb=100, bottom_rkb=500, diameter_in="12 1/4")
        assert d.diameter_in == 12.25

    def test_string_half_diameter(self):
        d = DrillingModel(top_rkb=325, bottom_rkb=720, diameter_in="17 1/2")
        assert d.diameter_in == 17.5

    def test_default_oh_perm(self):
        d = DrillingModel(top_rkb=100, bottom_rkb=500, diameter_in=12.25)
        assert d.oh_perm == 10000

    def test_invalid_diameter_raises(self):
        with pytest.raises(ValidationError):
            DrillingModel(top_rkb=100, bottom_rkb=500, diameter_in=[12])


# ── CasingCementModel ──────────────────────────────────────────────────


class TestCasingCementModel:
    def test_valid_casing(self):
        c = CasingCementModel(
            top_rkb=325,
            bottom_rkb=709,
            diameter_in="13 3/8",
            toc_rkb=325,
            boc_rkb=709,
            shoe=True,
        )
        assert c.diameter_in == 13.375
        assert c.shoe is True
        assert c.cb_perm == 5  # default

    def test_custom_cb_perm(self):
        c = CasingCementModel(
            top_rkb=325,
            bottom_rkb=709,
            diameter_in=13.375,
            toc_rkb=325,
            boc_rkb=709,
            shoe=True,
            cb_perm=0.1,
        )
        assert c.cb_perm == 0.1


# ── BarrierModel ────────────────────────────────────────────────────────


class TestBarrierModel:
    def test_valid_barrier(self):
        b = BarrierModel(
            barrier_name="cplug1",
            barrier_type="cplug",
            top_rkb=400,
            bottom_rkb=550,
        )
        assert b.barrier_name == "cplug1"
        assert b.barrier_perm == 0.5  # default

    def test_custom_perm(self):
        b = BarrierModel(
            barrier_name="cplug1",
            barrier_type="mech_plug",
            top_rkb=1000,
            bottom_rkb=1100,
            barrier_perm=5.0,
        )
        assert b.barrier_perm == 5.0


# ── GeologyModel ───────────────────────────────────────────────────────


class TestGeologyModel:
    def test_valid_geology(self):
        g = GeologyModel(top_rkb=325, geol_unit="OVERBURDEN", reservoir_flag=False)
        assert g.geol_unit == "OVERBURDEN"
        assert g.reservoir_flag is False

    def test_reservoir_flag_true(self):
        g = GeologyModel(top_rkb=1250, geol_unit="RESERVOIR", reservoir_flag=True)
        assert g.reservoir_flag is True


# ── ReservoirPressureModel ──────────────────────────────────────────────


class TestReservoirPressureModel:
    def test_with_string_pressures(self):
        rp = ReservoirPressureModel(depth_msl=1250, RP2="+ 25", RP3="- 25")
        assert rp.depth_msl == 1250
        assert rp.RP2 == "+ 25"
        assert rp.RP3 == "- 25"

    def test_optional_fields(self):
        rp = ReservoirPressureModel(depth_msl=2238)
        assert rp.RP1 is None
        assert rp.RP2 is None
        assert rp.RP3 is None


# ── MetaDataModel ───────────────────────────────────────────────────────


class TestMetaDataModel:
    def test_defaults(self):
        m = MetaDataModel()
        assert m.namespace == "screen"
        assert m.name is None
        assert m.author is None

    def test_with_values(self):
        m = MetaDataModel(namespace="prod", name="wildcat", author="equinor")
        assert m.namespace == "prod"


# ── WellSpec ────────────────────────────────────────────────────────────


class TestWellSpec:
    def test_minimal_spec(self):
        header = WellHeaderModel(
            well_name="TestWell",
            well_rkb=25,
            sf_depth_msl=300,
            well_td_rkb=1500,
            sf_temp=4,
            geo_tgrad=40,
        )
        spec = WellSpec(well_header=header, co2_datum=1300)
        assert spec.well_header.well_name == "TestWell"
        assert spec.drilling is None
        assert spec.barriers is None

    def test_with_drilling(self):
        header = WellHeaderModel(
            well_name="TestWell",
            well_rkb=25,
            sf_depth_msl=300,
            well_td_rkb=1500,
            sf_temp=4,
            geo_tgrad=40,
        )
        drilling = [
            DrillingModel(top_rkb=325, bottom_rkb=720, diameter_in="17 1/2"),
            DrillingModel(top_rkb=720, bottom_rkb=1185, diameter_in="12 1/4"),
        ]
        spec = WellSpec(well_header=header, drilling=drilling, co2_datum=1300)
        assert len(spec.drilling) == 2
        assert spec.drilling[0].diameter_in == 17.5


# ── WellModel (full YAML structure) ────────────────────────────────────


class TestWellModel:
    def test_full_model(self):
        model = WellModel(
            apiVersion="well/v0.1",
            kind="Well",
            metadata=MetaDataModel(namespace="screen", name="test", author="tester"),
            spec=WellPressureSpec(
                well_header=WellHeaderModel(
                    well_name="TestWell",
                    well_rkb=25,
                    sf_depth_msl=300,
                    well_td_rkb=1500,
                    sf_temp=4,
                    geo_tgrad=40,
                ),
                co2_datum=1300,
            ),
        )
        assert model.kind == "Well"
        assert model.spec.well_header.well_name == "TestWell"
