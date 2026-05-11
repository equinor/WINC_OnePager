"""Tests for Pressure, PressureScenarioManager, PressureScenario, and barrier_pressure modules."""

import json
import pathlib

import numpy as np
import pandas as pd
import pytest

from src.WellClass.libs.utils.yaml_parser import yaml_parser
from src.WellClass.libs.well_class.well_class import Well
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
            temps = below_sf["temperature"].to_numpy()
            assert temps[-1] > temps[0]

    def test_hydrostatic_pressure_increases_with_depth(self, pressure_instance):
        curves = pressure_instance.init_curves
        below_zero = curves[curves["depth"] > 0]
        pressures = below_zero["hydrostatic_pressure"].to_numpy()
        assert pressures[-1] > pressures[0]

    def test_shmin_increases_with_depth(self, pressure_instance):
        curves = pressure_instance.init_curves
        below_sf = curves[curves["depth"] > pressure_instance.sf_depth_msl]
        if len(below_sf) > 1:
            shmin = below_sf["min_horizontal_stress"].to_numpy()
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


# ── Pressure.add_scenario ──────────────────────────────────────────────


class TestPressureAddScenario:
    """Tests for add_scenario() with various parameter combinations, derived from notebook usage."""

    @pytest.fixture(scope="class")
    def pressure_base(self, simple_well_config):
        header = simple_well_config["well_header"]
        return Pressure(
            sf_depth_msl=header["sf_depth_msl"],
            well_td_rkb=header["well_td_rkb"],
            well_rkb=header["well_rkb"],
            sf_temp=header["sf_temp"],
            geo_tgrad=header["geo_tgrad"],
            fluid_type="pure_co2",
            pvt_path=PVT_PATH,
            z_fluid_contact=900.0,
            default_hs_scenario=False,
        )

    def test_add_scenario_p_delta_zero(self, pressure_base):
        pressure_base.add_scenario(
            scenario_name="hydrostatic_test",
            from_resrvr=True,
            p_delta=0,
        )
        sc = pressure_base.scenario_manager.scenarios["hydrostatic_test"]
        assert "fluid_pressure" in sc.init_curves.columns
        fp = sc.init_curves["fluid_pressure"]
        assert fp.notna().sum() > 0

    def test_add_scenario_p_delta_positive(self, pressure_base):
        pressure_base.add_scenario(
            scenario_name="overpressure",
            from_resrvr=True,
            p_delta=10.0,
        )
        sc = pressure_base.scenario_manager.scenarios["overpressure"]
        assert sc.p_delta == 10.0
        assert "fluid_pressure" in sc.init_curves.columns

    def test_add_scenario_p_delta_negative(self, pressure_base):
        pressure_base.add_scenario(
            scenario_name="underpressure",
            from_resrvr=True,
            p_delta=-5.0,
        )
        sc = pressure_base.scenario_manager.scenarios["underpressure"]
        assert sc.p_delta == -5.0

    def test_add_scenario_with_z_resrv_and_p_resrv(self, pressure_base):
        pressure_base.add_scenario(
            scenario_name="reservoir_point",
            from_resrvr=True,
            z_resrv=1000.0,
            p_resrv=105.0,
        )
        sc = pressure_base.scenario_manager.scenarios["reservoir_point"]
        assert sc.z_resrv == 1000.0
        assert sc.p_resrv == 105.0

    def test_add_scenario_with_specific_gravity(self, pressure_base):
        pressure_base.add_scenario(
            scenario_name="sg_scenario",
            from_resrvr=True,
            p_delta=5.0,
            specific_gravity=0.8,
        )
        sc = pressure_base.scenario_manager.scenarios["sg_scenario"]
        assert sc.specific_gravity == 0.8
        assert sc.fluid_interpolator is None

    def test_add_scenario_with_different_fluid_type(self, pressure_base):
        pressure_base.add_scenario(
            scenario_name="mixture1_scenario",
            from_resrvr=True,
            p_delta=0,
            fluid_type="mixture1",
        )
        sc = pressure_base.scenario_manager.scenarios["mixture1_scenario"]
        assert sc.fluid_type == "mixture1"
        assert sc.fluid_interpolator is not None

    def test_multiple_scenarios_coexist(self, pressure_base):
        names = pressure_base.scenario_manager.scenarios.keys()
        assert len(names) >= 3

    def test_duplicate_scenario_name_raises(self, pressure_base):
        with pytest.raises(ValueError, match="already exists"):
            pressure_base.add_scenario(
                scenario_name="hydrostatic_test",
                from_resrvr=True,
                p_delta=0,
            )


# ── Pressure init variants ─────────────────────────────────────────────


class TestPressureInitVariants:
    """Tests for Pressure instantiation with different parameter combinations seen in notebooks."""

    @pytest.fixture(scope="class")
    def header(self, simple_well_config):
        return simple_well_config["well_header"]

    def test_with_specific_gravity(self, header):
        p = Pressure(
            sf_depth_msl=header["sf_depth_msl"],
            well_td_rkb=header["well_td_rkb"],
            well_rkb=header["well_rkb"],
            sf_temp=header["sf_temp"],
            geo_tgrad=header["geo_tgrad"],
            fluid_type="pure_co2",
            pvt_path=PVT_PATH,
            specific_gravity=0.8,
        )
        assert p.fluid_interpolator is None
        assert p.brine_interpolator is not None

    def test_with_rho_brine(self, header):
        p = Pressure(
            sf_depth_msl=header["sf_depth_msl"],
            well_td_rkb=header["well_td_rkb"],
            well_rkb=header["well_rkb"],
            sf_temp=header["sf_temp"],
            geo_tgrad=header["geo_tgrad"],
            fluid_type="pure_co2",
            pvt_path=PVT_PATH,
            rho_brine=1025.0,
        )
        # With constant rho_brine, hydrostatic pressure should still be computed
        assert p.init_curves["hydrostatic_pressure"].notna().all()

    def test_with_shmin_gradient(self, header):
        custom_gradient = 0.2
        p = Pressure(
            sf_depth_msl=header["sf_depth_msl"],
            well_td_rkb=header["well_td_rkb"],
            well_rkb=header["well_rkb"],
            sf_temp=header["sf_temp"],
            geo_tgrad=header["geo_tgrad"],
            fluid_type="pure_co2",
            pvt_path=PVT_PATH,
            shmin_gradient=custom_gradient,
        )
        assert p.shmin_gradient == custom_gradient
        shmin = p.init_curves["min_horizontal_stress"]
        assert shmin.notna().sum() > 0

    def test_with_ip_shmin_data(self, header):
        shmin_data = np.array(
            [
                [0, 1.01325],
                [200, 25.0],
                [500, 60.0],
                [800, 95.0],
                [1200, 150.0],
            ]
        )
        p = Pressure(
            sf_depth_msl=header["sf_depth_msl"],
            well_td_rkb=header["well_td_rkb"],
            well_rkb=header["well_rkb"],
            sf_temp=header["sf_temp"],
            geo_tgrad=header["geo_tgrad"],
            fluid_type="pure_co2",
            pvt_path=PVT_PATH,
            ip_shmin_data=shmin_data,
        )
        shmin = p.init_curves["min_horizontal_stress"]
        assert shmin.notna().sum() > 0

    def test_default_hs_scenario_false(self, header):
        p = Pressure(
            sf_depth_msl=header["sf_depth_msl"],
            well_td_rkb=header["well_td_rkb"],
            well_rkb=header["well_rkb"],
            sf_temp=header["sf_temp"],
            geo_tgrad=header["geo_tgrad"],
            fluid_type="pure_co2",
            pvt_path=PVT_PATH,
            z_fluid_contact=900.0,
            default_hs_scenario=False,
        )
        assert len(p.scenario_manager.scenarios) == 0

    def test_default_hs_scenario_true(self, header):
        p = Pressure(
            sf_depth_msl=header["sf_depth_msl"],
            well_td_rkb=header["well_td_rkb"],
            well_rkb=header["well_rkb"],
            sf_temp=header["sf_temp"],
            geo_tgrad=header["geo_tgrad"],
            fluid_type="pure_co2",
            pvt_path=PVT_PATH,
            z_fluid_contact=900.0,
            default_hs_scenario=True,
        )
        assert "hydrostatic" in p.scenario_manager.scenarios

    def test_with_input_scenarios(self, header):
        scenarios = {"depth_msl": 900, "scenario_1": 5.0, "scenario_2": -3.0}
        p = Pressure(
            sf_depth_msl=header["sf_depth_msl"],
            well_td_rkb=header["well_td_rkb"],
            well_rkb=header["well_rkb"],
            sf_temp=header["sf_temp"],
            geo_tgrad=header["geo_tgrad"],
            fluid_type="pure_co2",
            pvt_path=PVT_PATH,
            z_fluid_contact=900.0,
            input_scenarios=scenarios,
        )
        assert "scenario_1" in p.scenario_manager.scenarios
        assert "scenario_2" in p.scenario_manager.scenarios
        assert len(p.scenario_manager.scenarios) == 2


# ── scenarios_summary ──────────────────────────────────────────────────


class TestScenariosSummary:
    def test_summary_has_expected_columns(self, simple_well_config):
        header = simple_well_config["well_header"]
        p = Pressure(
            sf_depth_msl=header["sf_depth_msl"],
            well_td_rkb=header["well_td_rkb"],
            well_rkb=header["well_rkb"],
            sf_temp=header["sf_temp"],
            geo_tgrad=header["geo_tgrad"],
            fluid_type="pure_co2",
            pvt_path=PVT_PATH,
            z_fluid_contact=900.0,
        )
        summary = p.scenarios_summary()
        assert isinstance(summary, pd.DataFrame)
        for col in ("name", "from_resrvr", "z_MSAD", "p_MSAD", "p_delta", "fluid_type"):
            assert col in summary.columns

    def test_summary_matches_manager(self, simple_well_config):
        header = simple_well_config["well_header"]
        p = Pressure(
            sf_depth_msl=header["sf_depth_msl"],
            well_td_rkb=header["well_td_rkb"],
            well_rkb=header["well_rkb"],
            sf_temp=header["sf_temp"],
            geo_tgrad=header["geo_tgrad"],
            fluid_type="pure_co2",
            pvt_path=PVT_PATH,
            z_fluid_contact=900.0,
        )
        s1 = p.scenarios_summary()
        s2 = p.scenario_manager.get_scenarios_summary()
        pd.testing.assert_frame_equal(s1, s2)


# ── Pressure.compute_barrier_leakage ───────────────────────────────────


class TestComputeBarrierLeakage:
    """Integration tests for compute_barrier_leakage, combining Pressure and Well."""

    @pytest.fixture(scope="class")
    def well_and_pressure(self):
        yaml_path = pathlib.Path("test_data/examples/wildcat/wildcat.yaml")
        well_model = yaml_parser(yaml_path)
        well_config = json.loads(well_model.spec.model_dump_json())
        header = well_config["well_header"]
        well = Well(
            header=well_config["well_header"],
            drilling=well_config["drilling"],
            casings=well_config["casing_cement"],
            geology=well_config["geology"],
            barriers=well_config["barriers"],
            barrier_perm=well_config.get("barrier_permeability"),
            co2_datum=well_config["co2_datum"],
        )

        p = Pressure(
            sf_depth_msl=header["sf_depth_msl"],
            well_td_rkb=header["well_td_rkb"],
            well_rkb=header["well_rkb"],
            sf_temp=header["sf_temp"],
            geo_tgrad=header["geo_tgrad"],
            fluid_type="pure_co2",
            pvt_path=PVT_PATH,
            z_fluid_contact=float(well_config["co2_datum"]),
        )
        return well, p

    def test_returns_dataframe(self, well_and_pressure):
        well, pressure = well_and_pressure
        barrier_name = next(iter(well.barriers_names.keys()))
        result = pressure.compute_barrier_leakage(well, barrier_name)
        assert isinstance(result, pd.DataFrame)

    def test_leakage_has_scenario_rows(self, well_and_pressure):
        well, pressure = well_and_pressure
        barrier_name = next(iter(well.barriers_names.keys()))
        result = pressure.compute_barrier_leakage(well, barrier_name)
        assert len(result) == len(pressure.scenario_manager.scenarios)

    def test_leakage_values_finite(self, well_and_pressure):
        well, pressure = well_and_pressure
        barrier_name = next(iter(well.barriers_names.keys()))
        result = pressure.compute_barrier_leakage(well, barrier_name)
        # Brine pressure above barrier should always be defined
        vals = pd.to_numeric(result["p_brine_above_barrier"])
        assert vals.notna().all()
        assert np.all(np.isfinite(vals))
