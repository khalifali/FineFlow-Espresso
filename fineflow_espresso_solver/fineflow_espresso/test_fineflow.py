"""Verification tests for the FineFlow-Espresso demonstrator."""

from __future__ import annotations

import unittest
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from fineflow_espresso import ModelParameters, _parameters_from_json, simulate


class FineFlowTests(unittest.TestCase):
    def short_case(self, **changes) -> ModelParameters:
        base = ModelParameters(
            control_mode="pressure",
            n_cells=24,
            duration_s=5.0,
            output_interval_s=0.25,
            max_time_step_s=0.05,
        )
        return replace(base, **changes)

    def assert_mass_conserved(self, result, tolerance=2e-11):
        initial = max(result.initial_fine_mass_kg, 1e-30)
        relative_error = np.max(np.abs(result.mass_balance_error_kg)) / initial
        self.assertLess(relative_error, tolerance)

    def test_pure_washthrough_conserves_mass(self):
        p = self.short_case(
            available_fines_initial_kg_m3=0.0,
            mobile_concentration_initial_kg_m3_liquid=2.0,
            release_rate_s=0.0,
            deposition_coefficient_m_inv=0.0,
            detachment_rate_s=0.0,
            axial_dispersion_m2_s=0.0,
        )
        result = simulate(p)
        self.assert_mass_conserved(result)
        self.assertGreater(result.cumulative_outlet_mass_kg[-1], 0.0)
        self.assertTrue(np.allclose(result.deposited_kg_m3, 0.0))

    def test_release_transfers_available_mass_without_creation(self):
        p = self.short_case(
            pressure_drop_pa=0.0,
            release_flow_exponent=0.0,
            deposition_coefficient_m_inv=0.0,
            detachment_rate_s=0.0,
            axial_dispersion_m2_s=0.0,
        )
        result = simulate(p)
        self.assert_mass_conserved(result)
        self.assertLess(
            np.mean(result.available_kg_m3[-1]),
            np.mean(result.available_kg_m3[0]),
        )
        self.assertGreater(np.mean(result.mobile_kg_m3_bulk[-1]), 0.0)

    def test_deposition_reduces_permeability_and_pressure_controlled_flow(self):
        p = self.short_case(
            duration_s=12.0,
            output_interval_s=0.5,
            release_rate_s=0.12,
            deposition_coefficient_m_inv=650.0,
            detachment_rate_s=0.0,
            axial_dispersion_m2_s=0.0,
        )
        result = simulate(p)
        self.assert_mass_conserved(result)
        self.assertGreater(np.max(result.deposited_kg_m3[-1]), 0.0)
        self.assertLess(
            np.min(result.permeability_m2[-1]),
            p.permeability_initial_m2,
        )
        self.assertLess(result.flow_rate_m3_s[-1], result.flow_rate_m3_s[0])

    def test_detachment_requires_threshold_exceedance(self):
        common = dict(
            available_fines_initial_kg_m3=0.0,
            deposited_fines_initial_kg_m3=8.0,
            mobile_concentration_initial_kg_m3_liquid=0.0,
            release_rate_s=0.0,
            deposition_coefficient_m_inv=0.0,
            axial_dispersion_m2_s=0.0,
            detachment_rate_s=0.8,
        )
        stable = simulate(
            self.short_case(critical_pressure_gradient_pa_m=1e15, **common)
        )
        eroding = simulate(
            self.short_case(critical_pressure_gradient_pa_m=1e5, **common)
        )
        self.assert_mass_conserved(stable)
        self.assert_mass_conserved(eroding)
        self.assertTrue(
            np.allclose(stable.deposited_kg_m3[-1], stable.deposited_kg_m3[0])
        )
        self.assertLess(
            np.mean(eroding.deposited_kg_m3[-1]),
            np.mean(eroding.deposited_kg_m3[0]),
        )

    def test_flow_control_predicts_pressure_increase(self):
        p = self.short_case(
            control_mode="flow",
            flow_rate_m3_s=1.0e-6,
            duration_s=8.0,
            release_rate_s=0.1,
            deposition_coefficient_m_inv=600.0,
            detachment_rate_s=0.0,
        )
        result = simulate(p)
        self.assert_mass_conserved(result)
        self.assertTrue(np.allclose(result.flow_rate_m3_s, p.flow_rate_m3_s))
        self.assertGreater(result.pressure_drop_pa[-1], result.pressure_drop_pa[0])

    def test_uniform_darcy_forchheimer_pressure_matches_exact_law(self):
        p = self.short_case(
            control_mode="flow",
            hydraulic_model="darcy_forchheimer",
            flow_rate_m3_s=5.0e-6,
            permeability_initial_m2=8.0e-13,
            forchheimer_coefficient_m_inv=1.2e5,
            basket_sieve_resistance_pa_s_m3=2.0e9,
            available_fines_initial_kg_m3=0.0,
            release_rate_s=0.0,
            deposition_coefficient_m_inv=0.0,
            detachment_rate_s=0.0,
        )
        result = simulate(p)
        q = p.flow_rate_m3_s / p.area_m2
        exact = p.length_m * (
            p.viscosity_pa_s * q / p.permeability_initial_m2
            + p.fluid_density_kg_m3
            * p.forchheimer_coefficient_m_inv
            * q**2
        ) + p.flow_rate_m3_s * p.basket_sieve_resistance_pa_s_m3
        self.assertTrue(np.allclose(result.pressure_drop_pa, exact, rtol=1e-12))

    def test_pressure_mode_inverts_darcy_forchheimer_law(self):
        p = self.short_case(
            control_mode="pressure",
            hydraulic_model="darcy_forchheimer",
            pressure_drop_pa=7.0e5,
            pressure_pulse_duration_s=0.0,
            permeability_initial_m2=8.0e-13,
            forchheimer_coefficient_m_inv=1.2e5,
            available_fines_initial_kg_m3=0.0,
            release_rate_s=0.0,
            deposition_coefficient_m_inv=0.0,
            detachment_rate_s=0.0,
        )
        result = simulate(p)
        self.assertTrue(np.allclose(result.pressure_drop_pa, p.pressure_drop_pa))

    def test_size_ratio_controls_successful_detachment_escape(self):
        common = dict(
            control_mode="flow",
            flow_rate_m3_s=1.0e-6,
            duration_s=1.0,
            available_fines_initial_kg_m3=0.0,
            deposited_fines_initial_kg_m3=8.0,
            release_rate_s=0.0,
            deposition_coefficient_m_inv=0.0,
            axial_dispersion_m2_s=0.0,
            detachment_rate_s=1.0,
            critical_pressure_gradient_pa_m=1.0e4,
            detachment_exponent=1.0,
            washout_escape_model="logistic",
            escape_ratio_midpoint=0.35,
            escape_ratio_steepness=25.0,
        )
        small = simulate(self.short_case(fine_particle_diameter_m=10e-6, **common))
        large = simulate(self.short_case(fine_particle_diameter_m=80e-6, **common))
        self.assertGreater(
            np.mean(large.deposited_kg_m3[-1]),
            np.mean(small.deposited_kg_m3[-1]),
        )
        self.assertGreater(
            np.mean(small.washout_escape_probability[0]),
            np.mean(large.washout_escape_probability[0]),
        )

    def test_pressure_pulse_erodes_blockage_and_increases_washout(self):
        p = self.short_case(
            duration_s=12.0,
            pressure_pulse_start_s=8.0,
            pressure_pulse_duration_s=2.0,
            pressure_pulse_pa=1.5e6,
            release_rate_s=0.15,
            deposition_coefficient_m_inv=700.0,
            detachment_rate_s=1.0,
            critical_pressure_gradient_pa_m=1.2e8,
        )
        result = simulate(p)
        self.assert_mass_conserved(result)
        before = np.argmin(np.abs(result.time_s - 7.75))
        pulse_start = np.argmin(np.abs(result.time_s - 8.0))
        pulse_end = np.argmin(np.abs(result.time_s - 10.0))
        pulse_mask = (result.time_s >= 8.0) & (result.time_s < 10.0)
        self.assertLess(
            result.deposited_kg_m3[pulse_end, -1],
            result.deposited_kg_m3[pulse_start, -1],
        )
        self.assertGreater(
            np.max(result.outlet_mass_rate_kg_s[pulse_mask]),
            result.outlet_mass_rate_kg_s[before],
        )

    def test_pi_controller_reaches_target_for_resistive_puck(self):
        p = self.short_case(
            control_mode="pi",
            duration_s=20.0,
            available_fines_initial_kg_m3=0.0,
            release_rate_s=0.0,
            deposition_coefficient_m_inv=0.0,
            detachment_rate_s=0.0,
        )
        result = simulate(p)
        self.assert_mass_conserved(result)
        self.assertAlmostEqual(
            result.pressure_drop_pa[-1] / 1e5,
            p.pressure_setpoint_pa / 1e5,
            delta=0.08,
        )
        self.assertLess(result.flow_rate_m3_s[-1], p.pump_flow_max_m3_s)

    def test_pi_flow_cap_leaves_coarse_puck_below_target(self):
        p = self.short_case(
            control_mode="pi",
            duration_s=15.0,
            permeability_initial_m2=1.0e-12,
            available_fines_initial_kg_m3=0.0,
            release_rate_s=0.0,
            deposition_coefficient_m_inv=0.0,
            detachment_rate_s=0.0,
        )
        result = simulate(p)
        self.assert_mass_conserved(result)
        self.assertLessEqual(
            np.max(result.flow_rate_m3_s), p.pump_flow_max_m3_s * (1 + 1e-12)
        )
        self.assertAlmostEqual(
            result.flow_rate_m3_s[-1], p.pump_flow_max_m3_s, delta=2e-9
        )
        self.assertLess(result.pressure_drop_pa[-1], 0.5 * p.pressure_setpoint_pa)
        self.assertTrue(result.controller_saturated[-1])

    def test_constant_pressure_ignores_pulse_and_pump_cap(self):
        p = self.short_case(
            control_mode="constant", hydraulic_model="darcy",
            permeability_initial_m2=1e-12,
            initial_permeability_profile_m2=tuple(np.linspace(1e-12, 2e-12, 24)),
            pressure_pulse_start_s=0.5, pressure_pulse_duration_s=2.0,
            basket_sieve_resistance_pa_s_m3=1e9,
            available_fines_initial_kg_m3=0.0,
            release_rate_s=0.0, deposition_coefficient_m_inv=0.0,
            detachment_rate_s=0.0,
        )
        result = simulate(p)
        resistance = p.viscosity_pa_s * p.dz_m * np.sum(
            1 / np.array(p.initial_permeability_profile_m2)) / p.area_m2
        expected_flow = p.pressure_drop_pa / (resistance + p.basket_sieve_resistance_pa_s_m3)
        np.testing.assert_allclose(result.pressure_drop_pa, p.pressure_drop_pa, rtol=1e-12)
        np.testing.assert_allclose(result.flow_rate_m3_s, expected_flow, rtol=1e-12)
        self.assertGreater(expected_flow, p.pump_flow_max_m3_s)
        self.assertFalse(np.any(result.controller_saturated))
        self.assert_mass_conserved(result)

    def test_constant_pressure_with_evolving_deposits(self):
        p = self.short_case(control_mode="constant", duration_s=10.0)
        result = simulate(p)
        np.testing.assert_allclose(result.pressure_drop_pa, p.pressure_drop_pa, rtol=1e-12)
        self.assertLess(result.flow_rate_m3_s[-1], result.flow_rate_m3_s[0])
        self.assert_mass_conserved(result)

    def test_named_case_config_and_prefixed_outputs(self):
        config = Path(__file__).with_name("case_config.json")
        case_name, parameters = _parameters_from_json(config)
        self.assertEqual(case_name, "fine_puck")
        self.assertEqual(parameters.control_mode, "constant")

        short_parameters = replace(
            parameters,
            duration_s=0.5,
            output_interval_s=0.5,
        )
        result = simulate(short_parameters)
        with TemporaryDirectory() as folder:
            result.save(folder, case_name="named test")
            names = {path.name for path in Path(folder).iterdir()}
            self.assertIn("named_test_timeseries.csv", names)
            self.assertIn("named_test_parameters.json", names)
            self.assertIn("named_test_dashboard.png", names)
            self.assertIn("named_test_hydraulic_response.png", names)
            self.assertIn("named_test_hydraulic_response.pdf", names)
            self.assertIn("named_test_fines_washout.png", names)
            self.assertIn("named_test_deposited_fines.png", names)
            self.assertIn("named_test_relative_permeability.png", names)

    def test_level_one_basket_splits_pressure_and_fines_mass(self):
        retention = 0.35
        sieve_resistance = 2.0e9
        p = self.short_case(
            control_mode="flow",
            flow_rate_m3_s=1.0e-6,
            basket_sieve_resistance_pa_s_m3=sieve_resistance,
            basket_fines_retention_fraction=retention,
            available_fines_initial_kg_m3=0.0,
            mobile_concentration_initial_kg_m3_liquid=2.0,
            release_rate_s=0.0,
            deposition_coefficient_m_inv=0.0,
            detachment_rate_s=0.0,
            axial_dispersion_m2_s=0.0,
        )
        result = simulate(p)
        self.assert_mass_conserved(result)
        self.assertTrue(
            np.allclose(
                result.basket_sieve_pressure_drop_pa,
                p.flow_rate_m3_s * sieve_resistance,
            )
        )
        self.assertTrue(
            np.allclose(
                result.pressure_drop_pa,
                result.puck_pressure_drop_pa
                + result.basket_sieve_pressure_drop_pa,
            )
        )
        self.assertAlmostEqual(
            result.cumulative_cup_mass_kg[-1],
            (1.0 - retention) * result.cumulative_outlet_mass_kg[-1],
            delta=1e-15,
        )
        self.assertAlmostEqual(
            result.cumulative_basket_retained_mass_kg[-1],
            retention * result.cumulative_outlet_mass_kg[-1],
            delta=1e-15,
        )

    def test_ct_profiles_and_pnm_permeability_table(self):
        n_cells = 8
        porosity_profile = tuple(np.linspace(0.34, 0.41, n_cells))
        permeability_profile = tuple(np.linspace(0.8e-14, 1.5e-14, n_cells))
        p = self.short_case(
            n_cells=n_cells,
            duration_s=0.25,
            output_interval_s=0.25,
            initial_porosity_profile=porosity_profile,
            initial_permeability_profile_m2=permeability_profile,
            deposited_fines_initial_kg_m3=5.0,
            available_fines_initial_kg_m3=0.0,
            release_rate_s=0.0,
            deposition_coefficient_m_inv=0.0,
            detachment_rate_s=0.0,
            permeability_model="pnm_table",
            pnm_deposit_kg_m3=(0.0, 5.0, 10.0),
            pnm_permeability_ratio=(1.0, 0.6, 0.2),
        )
        result = simulate(p)
        self.assertTrue(
            np.allclose(
                result.permeability_m2[0],
                np.asarray(permeability_profile) * 0.6,
            )
        )

    def test_pnm_hydraulic_diameter_and_escape_tables(self):
        p = self.short_case(
            n_cells=8,
            duration_s=0.25,
            output_interval_s=0.25,
            available_fines_initial_kg_m3=0.0,
            deposited_fines_initial_kg_m3=5.0,
            release_rate_s=0.0,
            deposition_coefficient_m_inv=0.0,
            detachment_rate_s=0.0,
            pnm_deposit_kg_m3=(0.0, 5.0, 10.0),
            pnm_hydraulic_diameter_ratio=(1.0, 0.8, 0.5),
            washout_escape_model="pnm_table",
            pnm_size_ratio=(0.0, 0.5, 1.0),
            pnm_escape_probability=(1.0, 0.4, 0.0),
        )
        result = simulate(p)
        epsilon0 = p.porosity_initial
        dh0 = (
            2.0
            / 3.0
            * p.particle_sphericity
            * p.coarse_effective_diameter_m
            * epsilon0
            / (1.0 - epsilon0)
        )
        expected_dh = 0.8 * dh0
        expected_ratio = p.fine_particle_diameter_m / expected_dh
        expected_escape = np.interp(
            expected_ratio, p.pnm_size_ratio, p.pnm_escape_probability
        )
        self.assertTrue(np.allclose(result.hydraulic_diameter_m[0], expected_dh))
        self.assertTrue(
            np.allclose(result.fine_to_hydraulic_diameter_ratio[0], expected_ratio)
        )
        self.assertTrue(
            np.allclose(result.washout_escape_probability[0], expected_escape)
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
