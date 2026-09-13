#!/usr/bin/env python3
"""FineFlow-Espresso: conservative 1D fines migration demonstrator.

The model resolves the espresso puck as axial finite-volume cells.  It tracks
fine mass in three states per unit bulk puck volume:

* ``available``: fines that may still be released from the coffee matrix,
* ``mobile``: fines suspended in the pore liquid,
* ``deposited``: fines retained in the pore space.

Local deposition reduces porosity and permeability.  A selectable Darcy or
Darcy-Forchheimer closure then updates the spatial pressure gradient and total
flow.  Hydraulic loading detaches deposits, while a size-ratio escape model
determines which detached fines can actually re-enter the mobile phase.

The default parameters are illustrative and are not fitted espresso data.
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ControlMode = Literal["constant", "pi", "pressure", "flow"]
PermeabilityModel = Literal["exponential", "kozeny_carman", "pnm_table"]
HydraulicModel = Literal["darcy", "darcy_forchheimer"]
WashoutEscapeModel = Literal["unity", "logistic", "pnm_table"]


@dataclass(frozen=True)
class ModelParameters:
    """Physical and effective model parameters in SI units."""

    # Geometry and fluid
    length_m: float = 0.020
    diameter_m: float = 0.058
    viscosity_pa_s: float = 1.0e-3
    fluid_density_kg_m3: float = 998.0
    n_cells: int = 40

    # Initial porous medium
    porosity_initial: float = 0.36
    permeability_initial_m2: float = 1.0e-14
    initial_porosity_profile: tuple[float, ...] = ()
    initial_permeability_profile_m2: tuple[float, ...] = ()
    deposit_density_kg_m3: float = 650.0
    minimum_porosity: float = 0.08
    minimum_permeability_ratio: float = 1.0e-5

    # Initial fines inventories (kg fines / m3 bulk puck)
    available_fines_initial_kg_m3: float = 28.0
    deposited_fines_initial_kg_m3: float = 0.0
    mobile_concentration_initial_kg_m3_liquid: float = 0.0
    deposit_capacity_kg_m3: float = 18.0

    # Release, transport, deposition, and detachment closures
    release_rate_s: float = 0.055
    release_flow_exponent: float = 0.70
    release_reference_velocity_m_s: float = 4.0e-4
    axial_dispersion_m2_s: float = 8.0e-8
    deposition_coefficient_m_inv: float = 420.0
    outlet_capture_multiplier: float = 5.0
    outlet_capture_decay_length_m: float = 0.0015
    detachment_rate_s: float = 0.80
    critical_pressure_gradient_pa_m: float = 1.35e8
    detachment_exponent: float = 1.5
    permeability_blocking_factor: float = 5.5
    permeability_model: PermeabilityModel = "exponential"
    pnm_deposit_kg_m3: tuple[float, ...] = ()
    pnm_permeability_ratio: tuple[float, ...] = ()

    # Hydraulic pressure-gradient closure.  For Darcy-Forchheimer,
    # G = mu*q/K + rho*beta_F*|q|*q.
    hydraulic_model: HydraulicModel = "darcy"
    forchheimer_coefficient_m_inv: float = 1.0e5
    initial_forchheimer_profile_m_inv: tuple[float, ...] = ()
    forchheimer_permeability_exponent: float = 0.0

    # Two-stage detachment and geometric escape/washout closure
    fine_particle_diameter_m: float = 25.0e-6
    coarse_effective_diameter_m: float = 300.0e-6
    particle_sphericity: float = 0.75
    initial_hydraulic_diameter_profile_m: tuple[float, ...] = ()
    pnm_hydraulic_diameter_ratio: tuple[float, ...] = ()
    washout_escape_model: WashoutEscapeModel = "logistic"
    escape_ratio_midpoint: float = 0.35
    escape_ratio_steepness: float = 18.0
    pnm_size_ratio: tuple[float, ...] = ()
    pnm_escape_probability: tuple[float, ...] = ()

    # Level-1 basket sieve: constant resistance and retention
    basket_sieve_resistance_pa_s_m3: float = 0.0
    basket_fines_retention_fraction: float = 0.0

    # Operating condition: constant never applies the legacy pressure pulse.
    # Pressure drop is across puck + basket, relative to the downstream outlet.
    control_mode: ControlMode = "constant"
    pressure_drop_pa: float = 9.0e5
    pressure_pulse_start_s: float = 50.0
    pressure_pulse_duration_s: float = 4.0
    pressure_pulse_pa: float = 1.50e6
    flow_rate_m3_s: float = 1.0e-6

    # Flow-actuated PI pressure controller
    pressure_setpoint_pa: float = 9.0e5
    pressure_setpoint_ramp_s: float = 2.0
    pump_flow_max_m3_s: float = 650.0e-6 / 60.0  # 650 mL/min = 10.83 mL/s
    pump_time_constant_s: float = 0.25
    pi_kp_ml_s_per_bar: float = 0.03
    pi_ki_ml_s_per_bar_s: float = 0.01
    controller_initial_flow_m3_s: float = 0.0

    # Integration and output
    duration_s: float = 60.0
    output_interval_s: float = 0.25
    max_time_step_s: float = 0.10
    cfl_safety: float = 0.42

    @property
    def area_m2(self) -> float:
        return np.pi * self.diameter_m**2 / 4.0

    @property
    def dz_m(self) -> float:
        return self.length_m / self.n_cells

    def validate(self) -> None:
        positive = {
            "length_m": self.length_m,
            "diameter_m": self.diameter_m,
            "viscosity_pa_s": self.viscosity_pa_s,
            "fluid_density_kg_m3": self.fluid_density_kg_m3,
            "n_cells": self.n_cells,
            "porosity_initial": self.porosity_initial,
            "permeability_initial_m2": self.permeability_initial_m2,
            "deposit_density_kg_m3": self.deposit_density_kg_m3,
            "deposit_capacity_kg_m3": self.deposit_capacity_kg_m3,
            "release_reference_velocity_m_s": self.release_reference_velocity_m_s,
            "outlet_capture_decay_length_m": self.outlet_capture_decay_length_m,
            "critical_pressure_gradient_pa_m": self.critical_pressure_gradient_pa_m,
            "fine_particle_diameter_m": self.fine_particle_diameter_m,
            "coarse_effective_diameter_m": self.coarse_effective_diameter_m,
            "particle_sphericity": self.particle_sphericity,
            "pressure_setpoint_pa": self.pressure_setpoint_pa,
            "pump_flow_max_m3_s": self.pump_flow_max_m3_s,
            "pump_time_constant_s": self.pump_time_constant_s,
            "duration_s": self.duration_s,
            "output_interval_s": self.output_interval_s,
            "max_time_step_s": self.max_time_step_s,
        }
        bad = [name for name, value in positive.items() if value <= 0]
        if bad:
            raise ValueError(f"Parameters must be positive: {', '.join(bad)}")
        if self.control_mode not in ("constant", "pi", "pressure", "flow"):
            raise ValueError("control_mode must be 'constant', 'pi', 'pressure', or 'flow'")
        if self.permeability_model not in (
            "exponential",
            "kozeny_carman",
            "pnm_table",
        ):
            raise ValueError(
                "permeability_model must be 'exponential', 'kozeny_carman', or 'pnm_table'"
            )
        if self.hydraulic_model not in ("darcy", "darcy_forchheimer"):
            raise ValueError("hydraulic_model must be 'darcy' or 'darcy_forchheimer'")
        if self.washout_escape_model not in ("unity", "logistic", "pnm_table"):
            raise ValueError(
                "washout_escape_model must be 'unity', 'logistic', or 'pnm_table'"
            )
        if not 0.0 <= self.basket_fines_retention_fraction <= 1.0:
            raise ValueError("basket_fines_retention_fraction must be between 0 and 1")
        for name, profile in (
            ("initial_porosity_profile", self.initial_porosity_profile),
            ("initial_permeability_profile_m2", self.initial_permeability_profile_m2),
            ("initial_forchheimer_profile_m_inv", self.initial_forchheimer_profile_m_inv),
            ("initial_hydraulic_diameter_profile_m", self.initial_hydraulic_diameter_profile_m),
        ):
            if profile and len(profile) != self.n_cells:
                raise ValueError(f"{name} must be empty or contain exactly n_cells values")
        if self.initial_porosity_profile and not all(
            self.minimum_porosity < value < 1.0
            for value in self.initial_porosity_profile
        ):
            raise ValueError(
                "initial_porosity_profile values must lie between minimum_porosity and 1"
            )
        if self.initial_permeability_profile_m2 and not all(
            value > 0.0 for value in self.initial_permeability_profile_m2
        ):
            raise ValueError("initial_permeability_profile_m2 values must be positive")
        if self.initial_forchheimer_profile_m_inv and not all(
            value >= 0.0 for value in self.initial_forchheimer_profile_m_inv
        ):
            raise ValueError("initial_forchheimer_profile_m_inv values must be non-negative")
        if self.initial_hydraulic_diameter_profile_m and not all(
            value > 0.0 for value in self.initial_hydraulic_diameter_profile_m
        ):
            raise ValueError("initial_hydraulic_diameter_profile_m values must be positive")
        if self.permeability_model == "pnm_table":
            if len(self.pnm_deposit_kg_m3) < 2 or len(self.pnm_deposit_kg_m3) != len(
                self.pnm_permeability_ratio
            ):
                raise ValueError(
                    "pnm_table requires equal deposit and permeability arrays with at least two points"
                )
            if np.any(np.diff(self.pnm_deposit_kg_m3) <= 0):
                raise ValueError("pnm_deposit_kg_m3 must be strictly increasing")
            if any(value <= 0 or value > 1 for value in self.pnm_permeability_ratio):
                raise ValueError("pnm_permeability_ratio values must lie in (0, 1]")
        if self.pnm_hydraulic_diameter_ratio:
            if len(self.pnm_deposit_kg_m3) < 2 or len(
                self.pnm_hydraulic_diameter_ratio
            ) != len(self.pnm_deposit_kg_m3):
                raise ValueError(
                    "pnm_hydraulic_diameter_ratio requires the same deposit grid as pnm_deposit_kg_m3"
                )
            if np.any(np.diff(self.pnm_deposit_kg_m3) <= 0):
                raise ValueError("pnm_deposit_kg_m3 must be strictly increasing")
            if any(value <= 0 or value > 1 for value in self.pnm_hydraulic_diameter_ratio):
                raise ValueError("pnm_hydraulic_diameter_ratio values must lie in (0, 1]")
        if self.washout_escape_model == "pnm_table":
            if len(self.pnm_size_ratio) < 2 or len(self.pnm_size_ratio) != len(
                self.pnm_escape_probability
            ):
                raise ValueError(
                    "pnm_table washout requires equal size-ratio and escape-probability arrays"
                )
            if np.any(np.diff(self.pnm_size_ratio) <= 0):
                raise ValueError("pnm_size_ratio must be strictly increasing")
            if any(value < 0 or value > 1 for value in self.pnm_escape_probability):
                raise ValueError("pnm_escape_probability values must lie in [0, 1]")
        if self.control_mode in ("constant", "pressure") and self.pressure_drop_pa < 0:
            raise ValueError("pressure_drop_pa must be non-negative")
        if self.control_mode == "flow" and self.flow_rate_m3_s < 0:
            raise ValueError("flow_rate_m3_s must be non-negative")
        if self.controller_initial_flow_m3_s > self.pump_flow_max_m3_s:
            raise ValueError("controller_initial_flow_m3_s cannot exceed pump_flow_max_m3_s")
        if not 0 < self.minimum_porosity < self.porosity_initial < 1:
            raise ValueError("Require 0 < minimum_porosity < porosity_initial < 1")
        if not 0 < self.minimum_permeability_ratio <= 1:
            raise ValueError("minimum_permeability_ratio must be in (0, 1]")
        lowest_initial_porosity = (
            min(self.initial_porosity_profile)
            if self.initial_porosity_profile
            else self.porosity_initial
        )
        maximum_safe_deposit = (
            lowest_initial_porosity - self.minimum_porosity
        ) * self.deposit_density_kg_m3
        if self.deposit_capacity_kg_m3 > maximum_safe_deposit:
            raise ValueError(
                "deposit_capacity_kg_m3 would reduce porosity below minimum_porosity"
            )
        nonnegative = {
            "available_fines_initial_kg_m3": self.available_fines_initial_kg_m3,
            "deposited_fines_initial_kg_m3": self.deposited_fines_initial_kg_m3,
            "mobile_concentration_initial_kg_m3_liquid": self.mobile_concentration_initial_kg_m3_liquid,
            "release_rate_s": self.release_rate_s,
            "release_flow_exponent": self.release_flow_exponent,
            "axial_dispersion_m2_s": self.axial_dispersion_m2_s,
            "deposition_coefficient_m_inv": self.deposition_coefficient_m_inv,
            "outlet_capture_multiplier": self.outlet_capture_multiplier,
            "detachment_rate_s": self.detachment_rate_s,
            "detachment_exponent": self.detachment_exponent,
            "permeability_blocking_factor": self.permeability_blocking_factor,
            "forchheimer_coefficient_m_inv": self.forchheimer_coefficient_m_inv,
            "forchheimer_permeability_exponent": self.forchheimer_permeability_exponent,
            "escape_ratio_midpoint": self.escape_ratio_midpoint,
            "escape_ratio_steepness": self.escape_ratio_steepness,
            "pressure_pulse_start_s": self.pressure_pulse_start_s,
            "pressure_pulse_duration_s": self.pressure_pulse_duration_s,
            "pressure_pulse_pa": self.pressure_pulse_pa,
            "pressure_setpoint_ramp_s": self.pressure_setpoint_ramp_s,
            "pi_kp_ml_s_per_bar": self.pi_kp_ml_s_per_bar,
            "pi_ki_ml_s_per_bar_s": self.pi_ki_ml_s_per_bar_s,
            "controller_initial_flow_m3_s": self.controller_initial_flow_m3_s,
            "basket_sieve_resistance_pa_s_m3": self.basket_sieve_resistance_pa_s_m3,
        }
        bad = [name for name, value in nonnegative.items() if value < 0]
        if bad:
            raise ValueError(f"Parameters must be non-negative: {', '.join(bad)}")
        if not 0.0 < self.particle_sphericity <= 1.0:
            raise ValueError("particle_sphericity must lie in (0, 1]")


@dataclass
class SimulationResult:
    """Time histories and full axial profiles returned by :func:`simulate`."""

    parameters: ModelParameters
    z_m: np.ndarray
    time_s: np.ndarray
    available_kg_m3: np.ndarray
    mobile_kg_m3_bulk: np.ndarray
    concentration_kg_m3_liquid: np.ndarray
    deposited_kg_m3: np.ndarray
    porosity: np.ndarray
    permeability_m2: np.ndarray
    forchheimer_coefficient_m_inv: np.ndarray
    pressure_gradient_pa_m: np.ndarray
    hydraulic_diameter_m: np.ndarray
    fine_to_hydraulic_diameter_ratio: np.ndarray
    washout_escape_probability: np.ndarray
    pressure_pa: np.ndarray
    flow_rate_m3_s: np.ndarray
    pump_flow_command_m3_s: np.ndarray
    pressure_drop_pa: np.ndarray
    puck_pressure_drop_pa: np.ndarray
    basket_sieve_pressure_drop_pa: np.ndarray
    pressure_setpoint_pa: np.ndarray
    controller_saturated: np.ndarray
    outlet_mass_rate_kg_s: np.ndarray
    cumulative_outlet_mass_kg: np.ndarray
    cup_mass_rate_kg_s: np.ndarray
    cumulative_cup_mass_kg: np.ndarray
    cumulative_basket_retained_mass_kg: np.ndarray
    mass_balance_error_kg: np.ndarray

    @property
    def initial_fine_mass_kg(self) -> float:
        p = self.parameters
        initial_density = (
            self.available_kg_m3[0]
            + self.mobile_kg_m3_bulk[0]
            + self.deposited_kg_m3[0]
        )
        return p.area_m2 * p.dz_m * np.sum(initial_density)

    def summary(self) -> dict[str, float | str]:
        p = self.parameters
        initial_k = _initial_permeability_profile(p)
        k_eff = p.length_m / np.sum(p.dz_m / self.permeability_m2[-1])
        blockage_index = int(np.argmin(self.permeability_m2[-1]))
        max_abs_error = float(np.max(np.abs(self.mass_balance_error_kg)))
        denominator = max(self.initial_fine_mass_kg, 1e-30)
        q = self.flow_rate_m3_s / p.area_m2
        inertial_drop = np.sum(
            p.fluid_density_kg_m3
            * self.forchheimer_coefficient_m_inv
            * (np.abs(q) * q)[:, None]
            * p.dz_m,
            axis=1,
        )
        inertial_fraction = np.divide(
            inertial_drop,
            self.puck_pressure_drop_pa,
            out=np.zeros_like(inertial_drop),
            where=np.abs(self.puck_pressure_drop_pa) > 0.0,
        )
        return {
            "control_mode": p.control_mode,
            "hydraulic_model": p.hydraulic_model,
            "permeability_model": p.permeability_model,
            "washout_escape_model": p.washout_escape_model,
            "duration_s": float(self.time_s[-1]),
            "initial_flow_rate_ml_s": float(self.flow_rate_m3_s[0] * 1e6),
            "final_flow_rate_ml_s": float(self.flow_rate_m3_s[-1] * 1e6),
            "maximum_flow_rate_ml_s": float(np.max(self.flow_rate_m3_s) * 1e6),
            "pump_flow_cap_ml_s": float(p.pump_flow_max_m3_s * 1e6),
            "initial_pressure_drop_bar": float(self.pressure_drop_pa[0] / 1e5),
            "final_pressure_drop_bar": float(self.pressure_drop_pa[-1] / 1e5),
            "maximum_pressure_drop_bar": float(np.max(self.pressure_drop_pa) / 1e5),
            "final_puck_pressure_drop_bar": float(self.puck_pressure_drop_pa[-1] / 1e5),
            "final_forchheimer_pressure_drop_bar": float(inertial_drop[-1] / 1e5),
            "maximum_forchheimer_fraction_of_puck_pressure": float(
                np.max(inertial_fraction)
            ),
            "final_basket_sieve_pressure_drop_bar": float(
                self.basket_sieve_pressure_drop_pa[-1] / 1e5
            ),
            "controller_saturated_fraction": float(np.mean(self.controller_saturated)),
            "final_effective_permeability_m2": float(k_eff),
            "minimum_local_permeability_ratio": float(
                np.min(self.permeability_m2[-1] / initial_k)
            ),
            "maximum_deposit_kg_m3": float(np.max(self.deposited_kg_m3[-1])),
            "maximum_pressure_gradient_MPa_m": float(
                np.max(self.pressure_gradient_pa_m) / 1e6
            ),
            "final_minimum_hydraulic_diameter_um": float(
                np.min(self.hydraulic_diameter_m[-1]) * 1e6
            ),
            "final_maximum_fine_to_hydraulic_diameter_ratio": float(
                np.max(self.fine_to_hydraulic_diameter_ratio[-1])
            ),
            "final_mean_washout_escape_probability": float(
                np.mean(self.washout_escape_probability[-1])
            ),
            "final_blockage_position_mm_from_inlet": float(
                self.z_m[blockage_index] * 1e3
            ),
            "cumulative_puck_outlet_fines_g": float(
                self.cumulative_outlet_mass_kg[-1] * 1e3
            ),
            "cumulative_basket_retained_fines_g": float(
                self.cumulative_basket_retained_mass_kg[-1] * 1e3
            ),
            "cumulative_cup_fines_g": float(self.cumulative_cup_mass_kg[-1] * 1e3),
            "maximum_relative_mass_balance_error": float(max_abs_error / denominator),
        }

    def timeseries_frame(self) -> pd.DataFrame:
        p = self.parameters
        min_k_ratio = np.min(
            self.permeability_m2 / _initial_permeability_profile(p), axis=1
        )
        max_gradient_mpa_m = np.max(self.pressure_gradient_pa_m, axis=1) / 1e6
        max_deposit = np.max(self.deposited_kg_m3, axis=1)
        blockage_index = np.argmin(self.permeability_m2, axis=1)
        initial_mass = max(self.initial_fine_mass_kg, 1e-30)
        q = self.flow_rate_m3_s / p.area_m2
        forchheimer_drop_pa = np.sum(
            p.fluid_density_kg_m3
            * self.forchheimer_coefficient_m_inv
            * (np.abs(q) * q)[:, None]
            * p.dz_m,
            axis=1,
        )
        viscous_drop_pa = self.puck_pressure_drop_pa - forchheimer_drop_pa
        forchheimer_fraction = np.divide(
            forchheimer_drop_pa,
            self.puck_pressure_drop_pa,
            out=np.zeros_like(forchheimer_drop_pa),
            where=np.abs(self.puck_pressure_drop_pa) > 0.0,
        )
        return pd.DataFrame(
            {
                "time_s": self.time_s,
                "flow_rate_ml_s": self.flow_rate_m3_s * 1e6,
                "pump_flow_command_ml_s": self.pump_flow_command_m3_s * 1e6,
                "pressure_drop_bar": self.pressure_drop_pa / 1e5,
                "puck_pressure_drop_bar": self.puck_pressure_drop_pa / 1e5,
                "viscous_puck_pressure_drop_bar": viscous_drop_pa / 1e5,
                "forchheimer_puck_pressure_drop_bar": forchheimer_drop_pa / 1e5,
                "forchheimer_fraction_of_puck_pressure": forchheimer_fraction,
                "basket_sieve_pressure_drop_bar": self.basket_sieve_pressure_drop_pa / 1e5,
                "pressure_setpoint_bar": self.pressure_setpoint_pa / 1e5,
                "controller_saturated": self.controller_saturated.astype(int),
                "puck_outlet_fines_mg_s": self.outlet_mass_rate_kg_s * 1e6,
                "cup_fines_mg_s": self.cup_mass_rate_kg_s * 1e6,
                "cumulative_puck_outlet_fines_g": self.cumulative_outlet_mass_kg * 1e3,
                "cumulative_basket_retained_fines_g": self.cumulative_basket_retained_mass_kg * 1e3,
                "cumulative_cup_fines_g": self.cumulative_cup_mass_kg * 1e3,
                "minimum_K_over_K0": min_k_ratio,
                "maximum_pressure_gradient_MPa_m": max_gradient_mpa_m,
                "maximum_deposit_kg_m3": max_deposit,
                "minimum_hydraulic_diameter_um": np.min(
                    self.hydraulic_diameter_m, axis=1
                )
                * 1e6,
                "maximum_fine_to_hydraulic_diameter_ratio": np.max(
                    self.fine_to_hydraulic_diameter_ratio, axis=1
                ),
                "mean_washout_escape_probability": np.mean(
                    self.washout_escape_probability, axis=1
                ),
                "blockage_position_mm_from_inlet": self.z_m[blockage_index] * 1e3,
                "relative_mass_balance_error": self.mass_balance_error_kg
                / initial_mass,
            }
        )

    def save(self, output_dir: Path | str, case_name: str = "default") -> None:
        """Save all outputs using a filesystem-safe case-name prefix."""

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        prefix = _safe_case_name(case_name)
        self.timeseries_frame().to_csv(
            output_dir / f"{prefix}_timeseries.csv", index=False
        )
        summary = {"case_name": case_name, **self.summary()}
        with (output_dir / f"{prefix}_summary.json").open(
            "w", encoding="utf-8"
        ) as stream:
            json.dump(summary, stream, indent=2)
        saved_config = {
            "case_name": case_name,
            "parameters": asdict(self.parameters),
        }
        with (output_dir / f"{prefix}_parameters.json").open(
            "w", encoding="utf-8"
        ) as stream:
            json.dump(saved_config, stream, indent=2)
        np.savez_compressed(
            output_dir / f"{prefix}_profiles.npz",
            z_m=self.z_m,
            time_s=self.time_s,
            available_kg_m3=self.available_kg_m3,
            mobile_kg_m3_bulk=self.mobile_kg_m3_bulk,
            concentration_kg_m3_liquid=self.concentration_kg_m3_liquid,
            deposited_kg_m3=self.deposited_kg_m3,
            porosity=self.porosity,
            permeability_m2=self.permeability_m2,
            forchheimer_coefficient_m_inv=self.forchheimer_coefficient_m_inv,
            pressure_gradient_pa_m=self.pressure_gradient_pa_m,
            hydraulic_diameter_m=self.hydraulic_diameter_m,
            fine_to_hydraulic_diameter_ratio=self.fine_to_hydraulic_diameter_ratio,
            washout_escape_probability=self.washout_escape_probability,
            pressure_pa=self.pressure_pa,
            total_pressure_drop_pa=self.pressure_drop_pa,
            puck_pressure_drop_pa=self.puck_pressure_drop_pa,
            basket_sieve_pressure_drop_pa=self.basket_sieve_pressure_drop_pa,
            puck_outlet_mass_rate_kg_s=self.outlet_mass_rate_kg_s,
            cup_mass_rate_kg_s=self.cup_mass_rate_kg_s,
            cumulative_puck_outlet_mass_kg=self.cumulative_outlet_mass_kg,
            cumulative_basket_retained_mass_kg=self.cumulative_basket_retained_mass_kg,
            cumulative_cup_mass_kg=self.cumulative_cup_mass_kg,
            pump_flow_command_m3_s=self.pump_flow_command_m3_s,
            pressure_setpoint_pa=self.pressure_setpoint_pa,
            controller_saturated=self.controller_saturated,
            superficial_velocity_m_s=np.repeat(
                (self.flow_rate_m3_s / self.parameters.area_m2)[:, None],
                self.parameters.n_cells,
                axis=1,
            ),
            pore_velocity_m_s=(self.flow_rate_m3_s / self.parameters.area_m2)[:, None]
            / self.porosity,
        )
        plot_result(self, output_dir / f"{prefix}_dashboard.png", case_name=case_name)
        plot_publication_figures(self, output_dir, case_name=case_name)


def _initial_porosity_profile(p: ModelParameters) -> np.ndarray:
    """Uniform default or a CT-derived axial porosity profile."""

    if p.initial_porosity_profile:
        return np.asarray(p.initial_porosity_profile, dtype=float)
    return np.full(p.n_cells, p.porosity_initial, dtype=float)


def _initial_permeability_profile(p: ModelParameters) -> np.ndarray:
    """Uniform default or a CT/PNM-derived axial permeability profile."""

    if p.initial_permeability_profile_m2:
        return np.asarray(p.initial_permeability_profile_m2, dtype=float)
    return np.full(p.n_cells, p.permeability_initial_m2, dtype=float)


def _initial_forchheimer_profile(p: ModelParameters) -> np.ndarray:
    """Uniform default or a CT/PNM/calibration-derived inertial profile."""

    if p.initial_forchheimer_profile_m_inv:
        return np.asarray(p.initial_forchheimer_profile_m_inv, dtype=float)
    return np.full(p.n_cells, p.forchheimer_coefficient_m_inv, dtype=float)


def _forchheimer_profile(
    p: ModelParameters, permeability: np.ndarray
) -> np.ndarray:
    """Return the current Forchheimer coefficient beta_F [1/m]."""

    if p.hydraulic_model == "darcy":
        return np.zeros(p.n_cells)
    k0 = _initial_permeability_profile(p)
    multiplier = (k0 / permeability) ** p.forchheimer_permeability_exponent
    return _initial_forchheimer_profile(p) * multiplier


def _initial_hydraulic_diameter_profile(p: ModelParameters) -> np.ndarray:
    """Return d_h,0 from CT/PNM input or the packed-bed surface relation."""

    if p.initial_hydraulic_diameter_profile_m:
        return np.asarray(p.initial_hydraulic_diameter_profile_m, dtype=float)
    epsilon0 = _initial_porosity_profile(p)
    return (
        (2.0 / 3.0)
        * p.particle_sphericity
        * p.coarse_effective_diameter_m
        * epsilon0
        / (1.0 - epsilon0)
    )


def _washout_geometry(
    p: ModelParameters, porosity: np.ndarray, deposited: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return hydraulic diameter, d_f/d_h, and successful escape probability."""

    dh0 = _initial_hydraulic_diameter_profile(p)
    if p.pnm_hydraulic_diameter_ratio:
        diameter_ratio = np.interp(
            deposited,
            np.asarray(p.pnm_deposit_kg_m3, dtype=float),
            np.asarray(p.pnm_hydraulic_diameter_ratio, dtype=float),
        )
    else:
        epsilon0 = _initial_porosity_profile(p)
        diameter_ratio = (
            porosity / (1.0 - porosity)
        ) / (epsilon0 / (1.0 - epsilon0))
    hydraulic_diameter = np.maximum(dh0 * diameter_ratio, 1.0e-12)
    size_ratio = p.fine_particle_diameter_m / hydraulic_diameter

    if p.washout_escape_model == "unity":
        escape_probability = np.ones_like(size_ratio)
    elif p.washout_escape_model == "logistic":
        argument = np.clip(
            p.escape_ratio_steepness * (size_ratio - p.escape_ratio_midpoint),
            -700.0,
            700.0,
        )
        escape_probability = 1.0 / (1.0 + np.exp(argument))
    else:
        escape_probability = np.interp(
            size_ratio,
            np.asarray(p.pnm_size_ratio, dtype=float),
            np.asarray(p.pnm_escape_probability, dtype=float),
        )
    return hydraulic_diameter, size_ratio, np.clip(escape_probability, 0.0, 1.0)


def _permeability_ratio(
    p: ModelParameters, deposited: np.ndarray, porosity: np.ndarray
) -> np.ndarray:
    """Evaluate the selected deposit-to-permeability closure."""

    if p.permeability_model == "exponential":
        ratio = np.exp(
            -p.permeability_blocking_factor
            * deposited
            / p.deposit_capacity_kg_m3
        )
    elif p.permeability_model == "kozeny_carman":
        epsilon0 = _initial_porosity_profile(p)
        ratio = (porosity / epsilon0) ** 3 * (
            (1.0 - epsilon0) / (1.0 - porosity)
        ) ** 2
    else:
        ratio = np.interp(
            deposited,
            np.asarray(p.pnm_deposit_kg_m3, dtype=float),
            np.asarray(p.pnm_permeability_ratio, dtype=float),
        )
    return np.clip(ratio, p.minimum_permeability_ratio, 1.0)


def _material_state(
    p: ModelParameters, deposited: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    porosity = _initial_porosity_profile(p) - deposited / p.deposit_density_kg_m3
    porosity = np.maximum(porosity, p.minimum_porosity)
    permeability = _initial_permeability_profile(p) * _permeability_ratio(
        p, deposited, porosity
    )
    return porosity, permeability


def _pressure_coefficients(
    p: ModelParameters, permeability: np.ndarray
) -> tuple[float, float, np.ndarray]:
    """Return total linear/quadratic coefficients for pressure as a function of q."""

    beta_f = _forchheimer_profile(p, permeability)
    linear = np.sum(p.viscosity_pa_s * p.dz_m / permeability)
    linear += p.area_m2 * p.basket_sieve_resistance_pa_s_m3
    quadratic = np.sum(p.fluid_density_kg_m3 * beta_f * p.dz_m)
    return float(linear), float(quadratic), beta_f


def _flow_rate_for_pressure(
    p: ModelParameters, permeability: np.ndarray, pressure_drop_pa: float
) -> float:
    """Invert the monotone Darcy-Forchheimer pressure law for Q >= 0."""

    if pressure_drop_pa <= 0.0:
        return 0.0
    linear, quadratic, _ = _pressure_coefficients(p, permeability)
    if quadratic <= 0.0:
        q = pressure_drop_pa / linear
    else:
        root = np.sqrt(linear**2 + 4.0 * quadratic * pressure_drop_pa)
        q = 2.0 * pressure_drop_pa / (linear + root)
    return float(q * p.area_m2)


def _hydraulics(
    p: ModelParameters,
    permeability: np.ndarray,
    time_s: float,
    controlled_flow_m3_s: float | None = None,
) -> tuple[float, float, np.ndarray, np.ndarray, float, float]:
    """Return velocity, total pressure, gradient, cell pressure, puck and sieve drops."""

    if p.control_mode == "pi":
        if controlled_flow_m3_s is None:
            raise ValueError("PI mode requires controlled_flow_m3_s")
        q = controlled_flow_m3_s / p.area_m2
    elif p.control_mode in ("constant", "pressure"):
        pulse_active = (
            p.control_mode == "pressure"
            and p.pressure_pulse_duration_s > 0
            and p.pressure_pulse_start_s <= time_s
            < p.pressure_pulse_start_s + p.pressure_pulse_duration_s
        )
        prescribed_pressure = (
            p.pressure_pulse_pa if pulse_active else p.pressure_drop_pa
        )
        flow_rate_m3_s = _flow_rate_for_pressure(
            p, permeability, prescribed_pressure
        )
        q = flow_rate_m3_s / p.area_m2
    else:
        q = p.flow_rate_m3_s / p.area_m2

    beta_f = _forchheimer_profile(p, permeability)
    gradient = (
        p.viscosity_pa_s * q / permeability
        + p.fluid_density_kg_m3 * beta_f * abs(q) * q
    )
    flow_rate_m3_s = q * p.area_m2
    sieve_pressure_drop = flow_rate_m3_s * p.basket_sieve_resistance_pa_s_m3
    puck_pressure_drop = float(np.sum(gradient * p.dz_m))
    pressure_drop = puck_pressure_drop + sieve_pressure_drop
    face_pressure = np.empty(p.n_cells + 1)
    face_pressure[0] = pressure_drop
    face_pressure[1:] = pressure_drop - np.cumsum(gradient * p.dz_m)
    face_pressure[-1] = sieve_pressure_drop  # pressure immediately above sieve
    cell_pressure = 0.5 * (face_pressure[:-1] + face_pressure[1:])
    return (
        q,
        pressure_drop,
        gradient,
        cell_pressure,
        puck_pressure_drop,
        sieve_pressure_drop,
    )


def _pressure_setpoint(p: ModelParameters, time_s: float) -> float:
    """Pressure target with an optional linear startup ramp."""

    if p.pressure_setpoint_ramp_s <= 0:
        return p.pressure_setpoint_pa
    return p.pressure_setpoint_pa * min(max(time_s / p.pressure_setpoint_ramp_s, 0.0), 1.0)


def _pi_controller_update(
    p: ModelParameters,
    dt: float,
    time_s: float,
    pressure_drop_pa: float,
    feedforward_flow_m3_s: float,
    actual_flow_m3_s: float,
    integral_error_bar_s: float,
) -> tuple[float, float, float, bool]:
    """Advance the flow actuator and a conditionally integrated PI controller."""

    setpoint_pa = _pressure_setpoint(p, time_s)
    error_bar = (setpoint_pa - pressure_drop_pa) / 1e5
    candidate_integral = integral_error_bar_s + error_bar * dt
    feedforward_ml_s = feedforward_flow_m3_s * 1e6
    raw_command_ml_s = (
        feedforward_ml_s
        + p.pi_kp_ml_s_per_bar * error_bar
        + p.pi_ki_ml_s_per_bar_s * candidate_integral
    )
    maximum_ml_s = p.pump_flow_max_m3_s * 1e6
    command_ml_s = float(np.clip(raw_command_ml_s, 0.0, maximum_ml_s))
    saturated_high = raw_command_ml_s > maximum_ml_s
    saturated_low = raw_command_ml_s < 0.0

    # Conditional integration prevents windup while the command is saturated.
    drives_out_of_saturation = (saturated_high and error_bar < 0) or (
        saturated_low and error_bar > 0
    )
    if (not saturated_high and not saturated_low) or drives_out_of_saturation:
        integral_error_bar_s = candidate_integral
        raw_command_ml_s = (
            feedforward_ml_s
            + p.pi_kp_ml_s_per_bar * error_bar
            + p.pi_ki_ml_s_per_bar_s * integral_error_bar_s
        )
        command_ml_s = float(np.clip(raw_command_ml_s, 0.0, maximum_ml_s))

    command_m3_s = command_ml_s * 1e-6
    actuator_fraction = -np.expm1(-dt / p.pump_time_constant_s)
    actual_flow_m3_s += (command_m3_s - actual_flow_m3_s) * actuator_fraction
    actual_flow_m3_s = float(
        np.clip(actual_flow_m3_s, 0.0, p.pump_flow_max_m3_s)
    )
    saturated = saturated_high or saturated_low
    return actual_flow_m3_s, command_m3_s, integral_error_bar_s, saturated


def _stable_time_step(
    p: ModelParameters, q: float, porosity: np.ndarray, remaining: float
) -> float:
    advection_rate = abs(q) / (np.min(porosity) * p.dz_m) if q else 0.0
    diffusion_rate = 2.0 * p.axial_dispersion_m2_s / p.dz_m**2
    transport_rate = advection_rate + diffusion_rate
    transport_limit = (
        p.cfl_safety / transport_rate if transport_rate > 0 else np.inf
    )
    controller_limit = 0.25 * p.pump_time_constant_s if p.control_mode == "pi" else np.inf
    return min(p.max_time_step_s, transport_limit, controller_limit, remaining)


def _local_reactions(
    p: ModelParameters,
    dt: float,
    q: float,
    gradient: np.ndarray,
    porosity: np.ndarray,
    available: np.ndarray,
    mobile: np.ndarray,
    deposited: np.ndarray,
    capture_multiplier: np.ndarray,
) -> None:
    """Apply positive, mass-conservative release/deposit/detach updates."""

    flow_ratio = abs(q) / p.release_reference_velocity_m_s
    release_multiplier = flow_ratio**p.release_flow_exponent
    released_fraction = -np.expm1(-p.release_rate_s * release_multiplier * dt)
    released = available * released_fraction
    available -= released
    mobile += released

    free_capacity_fraction = np.maximum(
        1.0 - deposited / p.deposit_capacity_kg_m3, 0.0
    )
    deposition_rate = (
        p.deposition_coefficient_m_inv
        * capture_multiplier
        * abs(q)
        * free_capacity_fraction
        / porosity
    )
    deposited_fraction = -np.expm1(-deposition_rate * dt)
    newly_deposited = np.minimum(
        mobile * deposited_fraction,
        np.maximum(p.deposit_capacity_kg_m3 - deposited, 0.0),
    )
    mobile -= newly_deposited
    deposited += newly_deposited

    # Hydraulic loading first detaches deposits.  The geometric escape
    # probability then determines which detached fines can actually enter the
    # mobile phase; the remainder is treated as immediately recaptured.
    _, _, escape_probability = _washout_geometry(p, porosity, deposited)
    excess_load = np.maximum(
        gradient / p.critical_pressure_gradient_pa_m - 1.0, 0.0
    )
    detachment_rate = p.detachment_rate_s * excess_load**p.detachment_exponent
    detached_fraction = -np.expm1(-detachment_rate * dt)
    successfully_remobilized = deposited * detached_fraction * escape_probability
    deposited -= successfully_remobilized
    mobile += successfully_remobilized


def _transport_mobile_fines(
    p: ModelParameters,
    dt: float,
    q: float,
    porosity: np.ndarray,
    mobile: np.ndarray,
) -> tuple[np.ndarray, float, float]:
    """Conservative upwind advection and centered axial dispersion.

    Clean water enters at z=0.  Both boundaries have zero diffusive fines flux.
    Only advective flux at z=L removes fines from the computational domain.
    """

    if q < 0:
        raise NotImplementedError("The demonstrator currently assumes inlet-to-outlet flow")

    concentration = mobile / porosity
    flux = np.zeros(p.n_cells + 1)
    flux[0] = 0.0  # clean inlet water
    flux[1:-1] = q * concentration[:-1]
    flux[-1] = q * concentration[-1]

    if p.axial_dispersion_m2_s > 0:
        porosity_face = 0.5 * (porosity[:-1] + porosity[1:])
        flux[1:-1] += (
            -porosity_face
            * p.axial_dispersion_m2_s
            * np.diff(concentration)
            / p.dz_m
        )

    updated = mobile - dt * np.diff(flux) / p.dz_m
    tolerance = 1e-12 * max(1.0, float(np.max(mobile)))
    if np.min(updated) < -tolerance:
        raise RuntimeError(
            "Negative mobile mass encountered; reduce max_time_step_s or cfl_safety"
        )
    updated = np.maximum(updated, 0.0)
    outlet_rate_kg_s = p.area_m2 * flux[-1]
    inlet_rate_kg_s = p.area_m2 * flux[0]
    return updated, outlet_rate_kg_s, inlet_rate_kg_s


def simulate(parameters: ModelParameters | None = None) -> SimulationResult:
    """Run a FineFlow simulation and return complete histories."""

    p = parameters or ModelParameters()
    p.validate()
    z = (np.arange(p.n_cells) + 0.5) * p.dz_m
    capture_multiplier = 1.0 + (p.outlet_capture_multiplier - 1.0) * np.exp(
        -(p.length_m - z) / p.outlet_capture_decay_length_m
    )

    available = np.full(
        p.n_cells, p.available_fines_initial_kg_m3, dtype=float
    )
    deposited = np.full(
        p.n_cells, p.deposited_fines_initial_kg_m3, dtype=float
    )
    porosity, permeability = _material_state(p, deposited)
    mobile = porosity * p.mobile_concentration_initial_kg_m3_liquid

    initial_total_mass = p.area_m2 * p.dz_m * np.sum(
        available + mobile + deposited
    )
    cumulative_outlet_mass = 0.0
    cumulative_cup_mass = 0.0
    cumulative_basket_retained_mass = 0.0
    cumulative_inlet_mass = 0.0
    controlled_flow_m3_s = (
        p.controller_initial_flow_m3_s if p.control_mode == "pi" else 0.0
    )
    pump_command_m3_s = controlled_flow_m3_s
    integral_error_bar_s = 0.0
    controller_is_saturated = False

    output_times = np.arange(
        0.0, p.duration_s + 0.5 * p.output_interval_s, p.output_interval_s
    )
    if output_times[-1] < p.duration_s:
        output_times = np.append(output_times, p.duration_s)
    else:
        output_times[-1] = p.duration_s
    n_output = len(output_times)

    shape = (n_output, p.n_cells)
    available_history = np.empty(shape)
    mobile_history = np.empty(shape)
    concentration_history = np.empty(shape)
    deposited_history = np.empty(shape)
    porosity_history = np.empty(shape)
    permeability_history = np.empty(shape)
    forchheimer_history = np.empty(shape)
    pressure_gradient_history = np.empty(shape)
    hydraulic_diameter_history = np.empty(shape)
    size_ratio_history = np.empty(shape)
    escape_probability_history = np.empty(shape)
    pressure_history = np.empty(shape)
    flow_history = np.empty(n_output)
    pump_command_history = np.empty(n_output)
    pressure_drop_history = np.empty(n_output)
    puck_pressure_drop_history = np.empty(n_output)
    basket_sieve_pressure_drop_history = np.empty(n_output)
    pressure_setpoint_history = np.empty(n_output)
    controller_saturated_history = np.empty(n_output, dtype=bool)
    outlet_rate_history = np.empty(n_output)
    cumulative_outlet_history = np.empty(n_output)
    cup_rate_history = np.empty(n_output)
    cumulative_cup_history = np.empty(n_output)
    cumulative_basket_retained_history = np.empty(n_output)
    mass_balance_history = np.empty(n_output)

    def record(index: int, outlet_rate: float) -> None:
        nonlocal porosity, permeability
        porosity, permeability = _material_state(p, deposited)
        (
            q,
            pressure_drop,
            gradient,
            pressure,
            puck_pressure_drop,
            sieve_pressure_drop,
        ) = _hydraulics(
            p,
            permeability,
            output_times[index],
            controlled_flow_m3_s if p.control_mode == "pi" else None,
        )
        available_history[index] = available
        mobile_history[index] = mobile
        concentration_history[index] = mobile / porosity
        deposited_history[index] = deposited
        porosity_history[index] = porosity
        permeability_history[index] = permeability
        forchheimer_history[index] = _forchheimer_profile(p, permeability)
        pressure_gradient_history[index] = gradient
        hydraulic_diameter, size_ratio, escape_probability = _washout_geometry(
            p, porosity, deposited
        )
        hydraulic_diameter_history[index] = hydraulic_diameter
        size_ratio_history[index] = size_ratio
        escape_probability_history[index] = escape_probability
        pressure_history[index] = pressure
        flow_history[index] = q * p.area_m2
        pump_command_history[index] = (
            pump_command_m3_s if p.control_mode == "pi" else q * p.area_m2
        )
        pressure_drop_history[index] = pressure_drop
        puck_pressure_drop_history[index] = puck_pressure_drop
        basket_sieve_pressure_drop_history[index] = sieve_pressure_drop
        pressure_setpoint_history[index] = (
            _pressure_setpoint(p, output_times[index])
            if p.control_mode == "pi"
            else (pressure_drop if p.control_mode in ("constant", "pressure") else np.nan)
        )
        controller_saturated_history[index] = (
            controller_is_saturated if p.control_mode == "pi" else False
        )
        outlet_rate_history[index] = outlet_rate
        cumulative_outlet_history[index] = cumulative_outlet_mass
        cup_rate_history[index] = (
            1.0 - p.basket_fines_retention_fraction
        ) * outlet_rate
        cumulative_cup_history[index] = cumulative_cup_mass
        cumulative_basket_retained_history[index] = cumulative_basket_retained_mass
        domain_mass = p.area_m2 * p.dz_m * np.sum(
            available + mobile + deposited
        )
        mass_balance_history[index] = (
            domain_mass
            + cumulative_cup_mass
            + cumulative_basket_retained_mass
            - cumulative_inlet_mass
            - initial_total_mass
        )

    time = 0.0
    most_recent_outlet_rate = 0.0
    record(0, most_recent_outlet_rate)

    for output_index in range(1, n_output):
        target_time = output_times[output_index]
        outlet_mass_since_output = 0.0
        interval_duration = target_time - time
        while time < target_time - 1e-14:
            porosity, permeability = _material_state(p, deposited)
            q, pressure_drop, gradient, _, _, _ = _hydraulics(
                p,
                permeability,
                time,
                controlled_flow_m3_s if p.control_mode == "pi" else None,
            )
            remaining = target_time - time
            if p.control_mode == "pressure" and p.pressure_pulse_duration_s > 0:
                pulse_events = (
                    p.pressure_pulse_start_s,
                    p.pressure_pulse_start_s + p.pressure_pulse_duration_s,
                )
                future_events = [event - time for event in pulse_events if event > time + 1e-14]
                if future_events:
                    remaining = min(remaining, min(future_events))
            dt = _stable_time_step(p, q, porosity, remaining)

            _local_reactions(
                p,
                dt,
                q,
                gradient,
                porosity,
                available,
                mobile,
                deposited,
                capture_multiplier,
            )
            porosity, _ = _material_state(p, deposited)
            mobile, outlet_rate, inlet_rate = _transport_mobile_fines(
                p, dt, q, porosity, mobile
            )
            outlet_increment = outlet_rate * dt
            inlet_increment = inlet_rate * dt
            cumulative_outlet_mass += outlet_increment
            retained_increment = (
                p.basket_fines_retention_fraction * outlet_increment
            )
            cup_increment = outlet_increment - retained_increment
            cumulative_basket_retained_mass += retained_increment
            cumulative_cup_mass += cup_increment
            cumulative_inlet_mass += inlet_increment
            outlet_mass_since_output += outlet_increment
            if p.control_mode == "pi":
                porosity, permeability = _material_state(p, deposited)
                _, pressure_drop, _, _, _, _ = _hydraulics(
                    p, permeability, time + dt, controlled_flow_m3_s
                )
                feedforward_flow_m3_s = _flow_rate_for_pressure(
                    p,
                    permeability,
                    _pressure_setpoint(p, time + dt),
                )
                (
                    controlled_flow_m3_s,
                    pump_command_m3_s,
                    integral_error_bar_s,
                    controller_is_saturated,
                ) = _pi_controller_update(
                    p,
                    dt,
                    time + dt,
                    pressure_drop,
                    feedforward_flow_m3_s,
                    controlled_flow_m3_s,
                    integral_error_bar_s,
                )
            time += dt

        most_recent_outlet_rate = outlet_mass_since_output / interval_duration
        record(output_index, most_recent_outlet_rate)

    return SimulationResult(
        parameters=p,
        z_m=z,
        time_s=output_times,
        available_kg_m3=available_history,
        mobile_kg_m3_bulk=mobile_history,
        concentration_kg_m3_liquid=concentration_history,
        deposited_kg_m3=deposited_history,
        porosity=porosity_history,
        permeability_m2=permeability_history,
        forchheimer_coefficient_m_inv=forchheimer_history,
        pressure_gradient_pa_m=pressure_gradient_history,
        hydraulic_diameter_m=hydraulic_diameter_history,
        fine_to_hydraulic_diameter_ratio=size_ratio_history,
        washout_escape_probability=escape_probability_history,
        pressure_pa=pressure_history,
        flow_rate_m3_s=flow_history,
        pump_flow_command_m3_s=pump_command_history,
        pressure_drop_pa=pressure_drop_history,
        puck_pressure_drop_pa=puck_pressure_drop_history,
        basket_sieve_pressure_drop_pa=basket_sieve_pressure_drop_history,
        pressure_setpoint_pa=pressure_setpoint_history,
        controller_saturated=controller_saturated_history,
        outlet_mass_rate_kg_s=outlet_rate_history,
        cumulative_outlet_mass_kg=cumulative_outlet_history,
        cup_mass_rate_kg_s=cup_rate_history,
        cumulative_cup_mass_kg=cumulative_cup_history,
        cumulative_basket_retained_mass_kg=cumulative_basket_retained_history,
        mass_balance_error_kg=mass_balance_history,
    )


def plot_result(
    result: SimulationResult,
    output_path: Path | str,
    case_name: str | None = None,
) -> None:
    """Create a four-panel diagnostic and prediction dashboard."""

    output_path = Path(output_path)
    time = result.time_s
    z_mm = result.z_m * 1e3
    p = result.parameters

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.4), constrained_layout=True)

    ax = axes[0, 0]
    flow_line = ax.plot(
        time,
        result.flow_rate_m3_s * 1e6,
        color="#005293",
        lw=2.1,
        label="Puck flow",
    )[0]
    legend_lines = [flow_line]
    if p.control_mode == "pi":
        command_line = ax.plot(
            time,
            result.pump_flow_command_m3_s * 1e6,
            color="#66B5D8",
            lw=1.3,
            ls=":",
            label="PI flow command",
        )[0]
        if np.any(result.controller_saturated):
            cap_line = ax.axhline(
                p.pump_flow_max_m3_s * 1e6,
                color="#005293",
                lw=0.9,
                ls="--",
                alpha=0.5,
                label="Pump-flow cap",
            )
            legend_lines.append(cap_line)
        legend_lines.append(command_line)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Flow rate [mL/s]", color="#005293")
    ax.tick_params(axis="y", labelcolor="#005293")
    ax.grid(alpha=0.22)
    pressure_axis = ax.twinx()
    pressure_line = pressure_axis.plot(
        time,
        result.pressure_drop_pa / 1e5,
        color="#E87722",
        lw=1.8,
        ls="--",
        label="Pressure drop",
    )[0]
    legend_lines.append(pressure_line)
    if p.control_mode == "pi":
        setpoint_line = pressure_axis.plot(
            time,
            result.pressure_setpoint_pa / 1e5,
            color="#B85D16",
            lw=1.0,
            ls=":",
            label="Pressure target",
        )[0]
        legend_lines.append(setpoint_line)
    pressure_axis.set_ylabel("Pressure drop [bar]", color="#B85D16")
    pressure_axis.tick_params(axis="y", labelcolor="#B85D16")
    ax.legend(handles=legend_lines, loc="best", fontsize=8, frameon=False)
    ax.set_title("Hydraulic response")

    ax = axes[0, 1]
    ax.plot(
        time,
        result.outlet_mass_rate_kg_s * 1e6,
        color="#E87722",
        lw=1.4,
        ls="--",
        label="Leaves puck",
    )
    ax.plot(
        time,
        result.cup_mass_rate_kg_s * 1e6,
        color="#B85D16",
        lw=2,
        label="Reaches cup",
    )
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Fines rate [mg/s]", color="#B85D16")
    ax.tick_params(axis="y", labelcolor="#B85D16")
    ax.grid(alpha=0.22)
    cumulative_axis = ax.twinx()
    cumulative_axis.plot(
        time, result.cumulative_cup_mass_kg * 1e3, color="#4C956C", lw=1.8
    )
    cumulative_axis.set_ylabel("Cumulative cup fines [g]", color="#36724E")
    cumulative_axis.tick_params(axis="y", labelcolor="#36724E")
    ax.set_title("Predicted fines washout")
    ax.legend(loc="best", fontsize=8, frameon=False)

    if p.control_mode == "pressure" and p.pressure_pulse_duration_s > 0:
        pulse_end = p.pressure_pulse_start_s + p.pressure_pulse_duration_s
        for pulse_axis in axes[0, :]:
            pulse_axis.axvspan(
                p.pressure_pulse_start_s,
                pulse_end,
                color="#E87722",
                alpha=0.10,
                lw=0,
            )

    extent = [time[0], time[-1], z_mm[-1] + 0.5 * p.dz_m * 1e3, z_mm[0] - 0.5 * p.dz_m * 1e3]
    ax = axes[1, 0]
    image = ax.imshow(
        result.deposited_kg_m3.T,
        aspect="auto",
        extent=extent,
        cmap="YlOrBr",
        interpolation="nearest",
    )
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Depth from inlet [mm]")
    ax.set_title("Deposited fines [kg/m$^3$ bulk]")
    fig.colorbar(image, ax=ax, pad=0.02)

    ax = axes[1, 1]
    permeability_ratio = np.maximum(
        result.permeability_m2 / _initial_permeability_profile(p), 1e-12
    )
    image = ax.imshow(
        np.log10(permeability_ratio).T,
        aspect="auto",
        extent=extent,
        cmap="viridis",
        interpolation="nearest",
        vmin=np.log10(p.minimum_permeability_ratio),
        vmax=0.0,
    )
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Depth from inlet [mm]")
    ax.set_title(r"Local permeability $\log_{10}(K/K_0)$")
    fig.colorbar(image, ax=ax, pad=0.02)

    case_label = f" | case: {case_name}" if case_name else ""
    fig.suptitle(
        f"FineFlow-Espresso{case_label} | illustrative, uncalibrated parameters",
        fontsize=14,
        fontweight="bold",
        color="#005293",
    )
    fig.savefig(output_path, dpi=180, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def _save_publication_figure(fig: plt.Figure, stem: Path) -> None:
    """Save one manuscript panel as high-resolution PNG and vector PDF."""

    png_path = stem.with_suffix(".png")
    pdf_path = stem.with_suffix(".pdf")
    temporary_png = stem.parent / f"{stem.name}.tmp.png"
    temporary_pdf = stem.parent / f"{stem.name}.tmp.pdf"
    fig.savefig(temporary_png, dpi=600, bbox_inches="tight", pad_inches=0.06)
    fig.savefig(temporary_pdf, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    temporary_png.replace(png_path)
    temporary_pdf.replace(pdf_path)


def plot_publication_figures(
    result: SimulationResult,
    output_dir: Path | str,
    case_name: str | None = None,
) -> None:
    """Save the four dashboard panels as standalone manuscript figures."""

    output_dir = Path(output_dir)
    prefix = _safe_case_name(case_name or "default")
    time = result.time_s
    z_mm = result.z_m * 1e3
    p = result.parameters
    extent = [
        time[0],
        time[-1],
        z_mm[-1] + 0.5 * p.dz_m * 1e3,
        z_mm[0] - 0.5 * p.dz_m * 1e3,
    ]

    fig, ax = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
    flow_line = ax.plot(
        time,
        result.flow_rate_m3_s * 1e6,
        color="#005293",
        lw=2.1,
        label="Puck flow",
    )[0]
    handles = [flow_line]
    if p.control_mode == "pi":
        command_line = ax.plot(
            time,
            result.pump_flow_command_m3_s * 1e6,
            color="#66B5D8",
            lw=1.3,
            ls=":",
            label="PI flow command",
        )[0]
        handles.append(command_line)
        if np.any(result.controller_saturated):
            handles.append(
                ax.axhline(
                    p.pump_flow_max_m3_s * 1e6,
                    color="#005293",
                    lw=0.9,
                    ls="--",
                    alpha=0.5,
                    label="Pump-flow cap",
                )
            )
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Flow rate [mL/s]", color="#005293")
    ax.tick_params(axis="y", labelcolor="#005293")
    ax.grid(alpha=0.22)
    pressure_axis = ax.twinx()
    handles.append(
        pressure_axis.plot(
            time,
            result.pressure_drop_pa / 1e5,
            color="#E87722",
            lw=1.8,
            ls="--",
            label="Pressure drop",
        )[0]
    )
    if p.control_mode == "pi":
        handles.append(
            pressure_axis.plot(
                time,
                result.pressure_setpoint_pa / 1e5,
                color="#B85D16",
                lw=1.0,
                ls=":",
                label="Pressure target",
            )[0]
        )
    pressure_axis.set_ylabel("Pressure drop [bar]", color="#B85D16")
    pressure_axis.tick_params(axis="y", labelcolor="#B85D16")
    ax.legend(handles=handles, loc="best", fontsize=8, frameon=False)
    ax.set_title("Hydraulic response")
    _save_publication_figure(fig, output_dir / f"{prefix}_hydraulic_response")

    fig, ax = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
    ax.plot(
        time,
        result.outlet_mass_rate_kg_s * 1e6,
        color="#E87722",
        lw=1.4,
        ls="--",
        label="Leaves puck",
    )
    ax.plot(
        time,
        result.cup_mass_rate_kg_s * 1e6,
        color="#B85D16",
        lw=2.0,
        label="Reaches cup",
    )
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Fines rate [mg/s]", color="#B85D16")
    ax.tick_params(axis="y", labelcolor="#B85D16")
    ax.grid(alpha=0.22)
    cumulative_axis = ax.twinx()
    cumulative_axis.plot(
        time,
        result.cumulative_cup_mass_kg * 1e3,
        color="#4C956C",
        lw=1.8,
    )
    cumulative_axis.set_ylabel("Cumulative cup fines [g]", color="#36724E")
    cumulative_axis.tick_params(axis="y", labelcolor="#36724E")
    ax.legend(loc="best", fontsize=8, frameon=False)
    ax.set_title("Predicted fines washout")
    _save_publication_figure(fig, output_dir / f"{prefix}_fines_washout")

    fig, ax = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
    image = ax.imshow(
        result.deposited_kg_m3.T,
        aspect="auto",
        extent=extent,
        cmap="YlOrBr",
        interpolation="nearest",
    )
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Depth from inlet [mm]")
    ax.set_title("Deposited fines [kg/m$^3$ bulk]")
    fig.colorbar(image, ax=ax, pad=0.02)
    _save_publication_figure(fig, output_dir / f"{prefix}_deposited_fines")

    permeability_ratio = np.maximum(
        result.permeability_m2 / _initial_permeability_profile(p), 1e-12
    )
    fig, ax = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
    image = ax.imshow(
        np.log10(permeability_ratio).T,
        aspect="auto",
        extent=extent,
        cmap="viridis",
        interpolation="nearest",
        vmin=np.log10(p.minimum_permeability_ratio),
        vmax=0.0,
    )
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Depth from inlet [mm]")
    ax.set_title(r"Local permeability $\log_{10}(K/K_0)$")
    fig.colorbar(image, ax=ax, pad=0.02)
    _save_publication_figure(fig, output_dir / f"{prefix}_relative_permeability")


def _safe_case_name(case_name: str) -> str:
    """Return a predictable filename-safe version of a user case name."""

    safe = re.sub(r"[^A-Za-z0-9_-]+", "_", case_name.strip()).strip("_")
    if not safe:
        raise ValueError("case_name must contain at least one letter or number")
    return safe


def _parameters_from_json(path: Path) -> tuple[str, ModelParameters]:
    """Read either a named multi-case config or a saved single-case config."""

    with path.open("r", encoding="utf-8") as stream:
        values = json.load(stream)

    if "selected_case" in values:
        case_name = str(values["selected_case"])
        cases = values.get("cases")
        if not isinstance(cases, dict) or case_name not in cases:
            raise ValueError(
                f"selected_case {case_name!r} is not present under 'cases'"
            )
        selected = cases[case_name]
        if not isinstance(selected, dict):
            raise ValueError(f"Case {case_name!r} must contain a JSON object")
        parameters = {k: v for k, v in selected.items() if not k.startswith("_")}
    elif "parameters" in values:
        case_name = str(values.get("case_name", path.stem))
        parameters = values["parameters"]
        if not isinstance(parameters, dict):
            raise ValueError("'parameters' must contain a JSON object")
    else:
        case_name = str(values.pop("case_name", path.stem))
        parameters = {k: v for k, v in values.items() if not k.startswith("_")}

    _safe_case_name(case_name)
    return case_name, replace(ModelParameters(), **parameters)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        help="JSON file containing ModelParameters fields; defaults are used otherwise.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results"),
        help="Directory for CSV, NPZ, JSON, PNG, and PDF outputs.",
    )
    parser.add_argument(
        "--case-name",
        help="Override the case name used in output filenames.",
    )
    parser.add_argument("--duration", type=float, help="Override duration [s].")
    parser.add_argument(
        "--pressure-bar", "--constant-pressure-bar", type=float,
        help="Hold this total puck + basket pressure drop [bar], without a pulse."
    )
    parser.add_argument(
        "--pi-pressure-bar",
        type=float,
        help="Use PI-controlled pressure with this target [bar].",
    )
    parser.add_argument(
        "--pump-cap-ml-s",
        type=float,
        help="Override the maximum pump flow used by PI control [mL/s].",
    )
    parser.add_argument(
        "--flow-ml-s", type=float, help="Use flow control at this volumetric flow [mL/s]."
    )
    parser.add_argument(
        "--hydraulic-model",
        choices=("darcy", "darcy_forchheimer"),
        help="Override the hydraulic pressure-gradient closure.",
    )
    parser.add_argument(
        "--forchheimer-coefficient-m-inv",
        type=float,
        help="Override the uniform Forchheimer coefficient beta_F [1/m].",
    )
    parser.add_argument(
        "--no-pulse",
        action="store_true",
        help="Disable the illustrative pressure pulse in pressure-control mode.",
    )
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    if args.config:
        case_name, parameters = _parameters_from_json(args.config)
    else:
        case_name, parameters = "default", ModelParameters()
    if args.case_name:
        case_name = args.case_name
        _safe_case_name(case_name)
    if args.duration is not None:
        parameters = replace(parameters, duration_s=args.duration)
    selected_modes = sum(
        value is not None
        for value in (args.pressure_bar, args.pi_pressure_bar, args.flow_ml_s)
    )
    if selected_modes > 1:
        raise SystemExit(
            "Choose only one of --pressure-bar, --pi-pressure-bar, or --flow-ml-s"
        )
    if args.pressure_bar is not None:
        parameters = replace(
            parameters,
            control_mode="constant",
            pressure_drop_pa=args.pressure_bar * 1e5,
        )
    if args.flow_ml_s is not None:
        parameters = replace(
            parameters,
            control_mode="flow",
            flow_rate_m3_s=args.flow_ml_s * 1e-6,
        )
    if args.pi_pressure_bar is not None:
        parameters = replace(
            parameters,
            control_mode="pi",
            pressure_setpoint_pa=args.pi_pressure_bar * 1e5,
        )
    if args.pump_cap_ml_s is not None:
        parameters = replace(
            parameters, pump_flow_max_m3_s=args.pump_cap_ml_s * 1e-6
        )
    if args.hydraulic_model is not None:
        parameters = replace(parameters, hydraulic_model=args.hydraulic_model)
    if args.forchheimer_coefficient_m_inv is not None:
        parameters = replace(
            parameters,
            forchheimer_coefficient_m_inv=args.forchheimer_coefficient_m_inv,
            initial_forchheimer_profile_m_inv=(),
        )
    if args.no_pulse:
        parameters = replace(parameters, pressure_pulse_duration_s=0.0)

    result = simulate(parameters)
    result.save(args.output_dir, case_name=case_name)
    print(json.dumps({"case_name": case_name, **result.summary()}, indent=2))
    print(f"\nSaved outputs to: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
