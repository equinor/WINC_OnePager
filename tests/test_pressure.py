"""Tests for Pressure, PressureScenarioManager, PressureScenario, and barrier_pressure modules."""

import json
import pathlib

import pandas as pd
import pytest

from src.WellClass.libs.utils.yaml_parser import yaml_parser
from src.WellClass.libs.well_pressure.barrier_pressure import leakage_proxy
from src.WellClass.libs.well_pressure.Pressure import Pressure
from src.WellClass.libs.well_pressure.PressureScenarioManager import PressureScenarioManager

PVT_PATH = pathlib.Path("src/WellClass/libs/pvt/pvt_constants")


@pytest.fixture(scope="module")
def simple_well_config():
    yaml_path = pathlib.Path("test_data/examples/simple_well/Simple_well.yaml")
    well_model = yaml_parser(yaml_path)
    return json.loads(well_model.spec.model_dump_json())


@pytest.fixture(scope="module")
def pressure_instance(simple_well_config):
    header = simple_well_config["well_header"]
    return Pressure(
        sf_depth_msl=header["sf_depth_msl"],
        well_td_rkb=header["well_td_rkb"],
        well_rkb=header["well_rkb"],
        sf_temp=header["sf_temp"],
        geo_tgrad=header["geo_tgrad"],
        fluid_type="pure_co2",
        pvt_path=PVT_PATH,
    )


# ── Pressure class ──────────────────────────────────────────────────────


class TestPressureInit:
    def test_pvt_data_loaded(self, pressure_instance):
        assert pressure_instance.pvt_data is not None
        assert "brine" in pressure_instance.pvt_data
        assert "pressure" in pressure_instance.pvt_data
        assert "temperature" in pressure_instance.pvt_data

    def test_interpolators_created(self, pressure_instance):
        assert pressure_instance.brine_interpolator is not None
        assert pressure_instance.fluid_interpolator is not None

    def test_init_curves_dataframe(self, pressure_instance):
        curves = pressure_instance.init_curves
        assert isinstance(curves, pd.DataFrame)
        for col in ["depth", "temperature", "hydrostatic_pressure", "min_horizontal_stress"]:
            assert col in curves.columns

    def test_depth_starts_at_or_above_zero(self, pressure_instance):
        depths = pressure_instance.init_curves["depth"]
        assert depths.iloc[0] <= 0

    def test_temperature_increases_with_depth(self, pressure_instance):
        curves = pressure_instance.init_curves
        # Below seafloor, temperature should increase
        below_sf = curves[curves["depth"] > pressure_instance.sf_depth_msl]
        if len(below_sf) > 1:
            temps = below_sf["temperature"].values
            assert temps[-1] > temps[0]

    def test_hydrostatic_pressure_increases_with_depth(self, pressure_instance):
        curves = pressure_instance.init_curves
        below_zero = curves[curves["depth"] > 0]
        pressures = below_zero["hydrostatic_pressure"].values
        assert pressures[-1] > pressures[0]

    def test_shmin_increases_with_depth(self, pressure_instance):
        curves = pressure_instance.init_curves
        below_sf = curves[curves["depth"] > pressure_instance.sf_depth_msl]
        if len(below_sf) > 1:
            shmin = below_sf["min_horizontal_stress"].values
            assert shmin[-1] > shmin[0]

    def test_scenario_manager_created(self, pressure_instance):
        assert isinstance(pressure_instance.scenario_manager, PressureScenarioManager)


# ── Pressure with scenarios ─────────────────────────────────────────────


class TestPressureWithScenarios:
    @pytest.fixture(scope="class")
    def pressure_with_scenarios(self, simple_well_config):
        header = simple_well_config["well_header"]
        rp = simple_well_config.get("reservoir_pressure", {})
        return Pressure(
            sf_depth_msl=header["sf_depth_msl"],
            well_td_rkb=header["well_td_rkb"],
            well_rkb=header["well_rkb"],
            sf_temp=header["sf_temp"],
            geo_tgrad=header["geo_tgrad"],
            fluid_type="pure_co2",
            pvt_path=PVT_PATH,
            z_fluid_contact=rp.get("depth_msl"),
            input_scenarios=rp,
        )

    def test_scenarios_created(self, pressure_with_scenarios):
        assert len(pressure_with_scenarios.scenario_manager.scenarios) > 0

    def test_scenarios_have_pressure_profiles(self, pressure_with_scenarios):
        for name, scenario in pressure_with_scenarios.scenario_manager.scenarios.items():
            assert "fluid_pressure" in scenario.init_curves.columns

    def test_scenarios_summary(self, pressure_with_scenarios):
        summary = pressure_with_scenarios.scenario_manager.get_scenarios_summary()
        assert isinstance(summary, pd.DataFrame)
        assert len(summary) > 0
        assert "name" in summary.columns


# ── PressureScenarioManager ────────────────────────────────────────────


class TestPressureScenarioManager:
    def test_create_scenario(self, pressure_instance):
        manager = PressureScenarioManager()
        scenario = manager.create_scenario(
            name="test_scenario",
            from_resrvr=True,
            init_curves=pressure_instance.init_curves,
            brine_interpolator=pressure_instance.brine_interpolator,
            fluid_interpolator=pressure_instance.fluid_interpolator,
            fluid_type="pure_co2",
            pvt_data=pressure_instance.pvt_data,
            z_resrv=1200.0,
            p_resrv=120.0,
            z_fluid_contact=1250.0,
        )
        assert scenario.name == "test_scenario"
        assert "test_scenario" in manager.scenarios

    def test_duplicate_name_raises(self, pressure_instance):
        manager = PressureScenarioManager()
        manager.create_scenario(
            name="dup",
            from_resrvr=True,
            init_curves=pressure_instance.init_curves,
            brine_interpolator=pressure_instance.brine_interpolator,
        )
        with pytest.raises(ValueError, match="already exists"):
            manager.create_scenario(
                name="dup",
                from_resrvr=True,
                init_curves=pressure_instance.init_curves,
                brine_interpolator=pressure_instance.brine_interpolator,
            )


# ── leakage_proxy ───────────────────────────────────────────────────────


class TestLeakageProxy:
    def test_zero_pressure_difference(self):
        """No pressure difference → no leakage (or clipped to zero)."""
        result = leakage_proxy(
            rho_fluid_below_barrier=800.0,
            rho_brine_below_barrier=1025.0,
            p_fluid_below_barrier=100.0,
            p_brine_above_barrier=100.0,
            permeability=0.01,
            barrier_props={"radius": 0.1, "height": 100.0},
        )
        assert result >= 0

    def test_higher_permeability_more_leakage(self):
        """Higher permeability should result in equal or more leakage."""
        common = {
            "rho_fluid_below_barrier": 600.0,
            "rho_brine_below_barrier": 1025.0,
            "p_fluid_below_barrier": 150.0,
            "p_brine_above_barrier": 100.0,
            "barrier_props": {"radius": 0.1, "height": 50.0},
        }
        low_perm = leakage_proxy(permeability=0.001, **common)
        high_perm = leakage_proxy(permeability=1.0, **common)
        assert high_perm >= low_perm

    def test_leakage_is_non_negative(self):
        """Leakage should never be negative."""
        result = leakage_proxy(
            rho_fluid_below_barrier=1025.0,
            rho_brine_below_barrier=1025.0,
            p_fluid_below_barrier=50.0,
            p_brine_above_barrier=200.0,
            permeability=0.01,
            barrier_props={"radius": 0.1, "height": 100.0},
        )
        assert result >= 0
