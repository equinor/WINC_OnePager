"""Tests for PVT module: data loading, interpolation, density, envelopes, and integration."""

import pathlib

import numpy as np
import pandas as pd
import pytest
import scipy.constants as const
from scipy.interpolate import RectBivariateSpline

from src.WellClass.libs.pvt.pvt import (
    _check_phase,
    _eval_envelopes_for_temps,
    _integrate_pressure,
    build_ivp_data,
    corr_rhobrine_LaliberteCopper,
    get_rho_from_pvt_data,
    load_envelopes,
    load_pvt_data,
    set_envelope_interpolator,
    solve_ivp_with_data,
)

PVT_PATH = pathlib.Path("src/WellClass/libs/pvt/pvt_constants")


# ── load_pvt_data ───────────────────────────────────────────────────────


class TestLoadPvtData:
    def test_brine_only(self):
        data = load_pvt_data(PVT_PATH, load_fluid=False)
        assert "pressure" in data
        assert "temperature" in data
        assert "brine" in data
        assert "rho" in data["brine"]

    def test_pure_co2(self):
        data = load_pvt_data(PVT_PATH, fluid_type="pure_co2", load_fluid=True)
        assert "pure_co2" in data
        assert "rho" in data["pure_co2"]
        assert "metadata" in data["pure_co2"]

    def test_mixture1(self):
        data = load_pvt_data(PVT_PATH, fluid_type="mixture1", load_fluid=True)
        assert "mixture1" in data
        assert "rho" in data["mixture1"]

    def test_pressure_temperature_are_1d_arrays(self):
        data = load_pvt_data(PVT_PATH, load_fluid=False)
        assert data["pressure"].ndim == 1
        assert data["temperature"].ndim == 1

    def test_rho_is_2d_matrix(self):
        data = load_pvt_data(PVT_PATH, fluid_type="pure_co2", load_fluid=True)
        assert data["brine"]["rho"].ndim == 2
        assert data["pure_co2"]["rho"].ndim == 2

    def test_invalid_fluid_type_raises(self):
        with pytest.raises(FileNotFoundError):
            load_pvt_data(PVT_PATH, fluid_type="nonexistent_fluid", load_fluid=True)

    def test_metadata_has_composition(self):
        data = load_pvt_data(PVT_PATH, fluid_type="pure_co2", load_fluid=True)
        assert "composition" in data["pure_co2"]["metadata"]


# ── load_envelopes / set_envelope_interpolator ──────────────────────────


class TestEnvelopes:
    @pytest.fixture(scope="class")
    def co2_envelopes(self):
        return load_envelopes(PVT_PATH / "pure_co2")

    def test_returns_bubble_and_dew(self, co2_envelopes):
        assert "bubble_point" in co2_envelopes
        assert "dew_point" in co2_envelopes

    def test_envelope_data_is_2d(self, co2_envelopes):
        for key in ("bubble_point", "dew_point"):
            arr = co2_envelopes[key]
            assert arr.ndim == 2
            assert arr.shape[1] == 2  # (T, P) columns

    def test_set_envelope_interpolator_callable(self, co2_envelopes):
        interp = set_envelope_interpolator(co2_envelopes["bubble_point"])
        assert callable(interp)

    def test_set_envelope_interpolator_none_input(self):
        result = set_envelope_interpolator(None)
        assert result is None

    def test_interpolator_returns_finite_value(self, co2_envelopes):
        interp = set_envelope_interpolator(co2_envelopes["bubble_point"])
        T_mid = co2_envelopes["bubble_point"][:, 0].mean()
        result = float(interp(T_mid))
        assert np.isfinite(result)


# ── get_rho_from_pvt_data ───────────────────────────────────────────────


class TestGetRhoFromPvtData:
    @pytest.fixture(scope="class")
    def brine_interpolator(self):
        data = load_pvt_data(PVT_PATH, load_fluid=False)
        return RectBivariateSpline(data["pressure"], data["temperature"], data["brine"]["rho"])

    def test_returns_positive_density(self, brine_interpolator):
        rho = get_rho_from_pvt_data(pressure=100.0, temperature=20.0, rho_interpolator=brine_interpolator)
        assert rho > 0

    def test_density_in_reasonable_range(self, brine_interpolator):
        rho = get_rho_from_pvt_data(pressure=200.0, temperature=50.0, rho_interpolator=brine_interpolator)
        assert 900 < rho < 1200  # brine density range kg/m³


# ── corr_rhobrine_LaliberteCopper ───────────────────────────────────────


class TestBrineDensityCorrection:
    def test_returns_array(self):
        temperature = np.array([10.0, 20.0, 30.0])
        pressure = np.array([100.0, 200.0, 300.0])
        rho_h2o = np.array([1000.0, 998.0, 995.0])
        result = corr_rhobrine_LaliberteCopper(salinity=3.5, temperature=temperature, pressure=pressure, rho_h2o=rho_h2o)
        assert result.shape == temperature.shape

    def test_brine_denser_than_freshwater(self):
        temperature = np.array([15.0, 25.0])
        pressure = np.array([100.0, 200.0])
        rho_h2o = np.array([999.0, 997.0])
        result = corr_rhobrine_LaliberteCopper(salinity=3.5, temperature=temperature, pressure=pressure, rho_h2o=rho_h2o)
        assert np.all(result > rho_h2o)

    def test_zero_salinity_equals_freshwater(self):
        temperature = np.array([20.0])
        pressure = np.array([100.0])
        rho_h2o = np.array([998.0])
        result = corr_rhobrine_LaliberteCopper(salinity=0.0, temperature=temperature, pressure=pressure, rho_h2o=rho_h2o)
        np.testing.assert_allclose(result, rho_h2o, rtol=1e-3)


# ── _check_phase ────────────────────────────────────────────────────────


class TestCheckPhase:
    def test_supercritical_above_tcrit(self):
        assert _check_phase(P=100, T=50, T_crit=31, bubble_interp=None, dew_interp=None) == "supercritical"

    def test_unknown_without_envelopes(self):
        assert _check_phase(P=100, T=20, T_crit=None, bubble_interp=None, dew_interp=None) == "unknown"


# ── _eval_envelopes_for_temps ───────────────────────────────────────────


class TestEvalEnvelopes:
    def test_returns_nan_without_interpolators(self):
        temps = np.array([10.0, 20.0, 30.0])
        p_dew, p_bub = _eval_envelopes_for_temps(temps, dew_interp=None, bub_interp=None)
        assert np.all(np.isnan(p_dew))
        assert np.all(np.isnan(p_bub))

    def test_invalidates_above_tcrit(self):
        envelopes = load_envelopes(PVT_PATH / "pure_co2")
        bub_interp = set_envelope_interpolator(envelopes["bubble_point"])
        dew_interp = set_envelope_interpolator(envelopes["dew_point"])

        temps = np.array([10.0, 20.0, 50.0])  # 50°C > T_crit for CO2 (31.1°C)
        p_dew, p_bub = _eval_envelopes_for_temps(temps, dew_interp, bub_interp, T_crit=31.1)
        assert np.isnan(p_dew[2])
        assert np.isnan(p_bub[2])


# ── _integrate_pressure ─────────────────────────────────────────────────


class TestIntegratePressure:
    @pytest.fixture(scope="class")
    def pvt_setup(self):
        data = load_pvt_data(PVT_PATH, fluid_type="pure_co2", load_fluid=True)
        brine_interp = RectBivariateSpline(data["pressure"], data["temperature"], data["brine"]["rho"])
        depth = np.linspace(0, 2000, 201)
        init_curves = pd.DataFrame(
            {
                "depth": depth,
                "temperature": 4 + 0.03 * depth,
                "hydrostatic_pressure": np.nan,
                "min_horizontal_stress": np.nan,
            }
        )
        return data, brine_interp, init_curves

    def test_pressure_monotonically_increases(self, pvt_setup):
        data, interp, curves = pvt_setup

        pressures, _ = _integrate_pressure(
            init_curves=curves,
            reference_depth=0,
            reference_pressure=const.atm / const.bar,
            pvt_data=data,
            fluid_key="brine",
            interpolator=interp,
        )
        finite = pressures[np.isfinite(pressures)]
        assert len(finite) > 10
        assert np.all(np.diff(finite) >= 0)

    def test_returns_warnings_list(self, pvt_setup):
        data, interp, curves = pvt_setup

        _, warnings = _integrate_pressure(
            init_curves=curves,
            reference_depth=0,
            reference_pressure=const.atm / const.bar,
            pvt_data=data,
            fluid_key="brine",
            interpolator=interp,
        )
        assert isinstance(warnings, list)


# ── build_ivp_data / solve_ivp_with_data ────────────────────────────────


class TestIvpHelpers:
    def test_build_ivp_data_returns_expected_keys(self):
        depth = np.linspace(0, 1000, 101)
        curves = pd.DataFrame({"depth": depth, "temperature": 4 + 0.03 * depth})
        pvt_data = load_pvt_data(PVT_PATH, fluid_type="pure_co2", load_fluid=True)
        result = build_ivp_data(500, 0, 50.0, curves, pvt_data["brine"])
        for key in ("interval_start", "interval_end", "initial_state", "index", "t_eval", "depth_array", "temperature_array"):
            assert key in result

    def test_solve_ivp_produces_solution(self):
        pvt_data = load_pvt_data(PVT_PATH, fluid_type="pure_co2", load_fluid=True)
        interp = RectBivariateSpline(pvt_data["pressure"], pvt_data["temperature"], pvt_data["brine"]["rho"])
        depth = np.linspace(0, 1000, 101)
        curves = pd.DataFrame({"depth": depth, "temperature": 4 + 0.03 * depth})
        ivp_data = build_ivp_data(0, 1000, 1.01325, curves, pvt_data["brine"])
        sol = solve_ivp_with_data(ivp_data, interp, "brine")
        assert sol.success
        assert len(sol.y[0]) > 0
