#!/usr/bin/env python3
"""Independent numerical-verification suite for FineFlow-Espresso.

This program complements the unit tests. It compares selected model components
with analytical solutions, checks invariants and boundedness, and performs
spatial- and temporal-refinement studies. It is verification of the numerical
implementation, not validation of the physical closures against experiments.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fineflow_espresso import ModelParameters, simulate


@dataclass
class Check:
    name: str
    passed: bool
    value: float
    limit: float
    units: str = ""
    explanation: str = ""


def inert_case(**changes) -> ModelParameters:
    """Return a case with all fines source terms disabled."""

    base = ModelParameters(
        control_mode="flow",
        hydraulic_model="darcy",
        flow_rate_m3_s=1.0e-6,
        duration_s=2.0,
        output_interval_s=0.1,
        max_time_step_s=0.02,
        available_fines_initial_kg_m3=0.0,
        deposited_fines_initial_kg_m3=0.0,
        mobile_concentration_initial_kg_m3_liquid=0.0,
        release_rate_s=0.0,
        deposition_coefficient_m_inv=0.0,
        detachment_rate_s=0.0,
        axial_dispersion_m2_s=0.0,
    )
    return replace(base, **changes)


def relative_error(value: float, reference: float, scale: float = 1e-30) -> float:
    return abs(value - reference) / max(abs(reference), scale)


def hydraulic_and_velocity_checks() -> tuple[list[Check], dict[str, np.ndarray]]:
    p = inert_case(duration_s=0.5, output_interval_s=0.1)
    result = simulate(p)
    exact_dp = p.viscosity_pa_s * p.length_m * (
        p.flow_rate_m3_s / p.area_m2
    ) / p.permeability_initial_m2
    darcy_error = float(
        np.max(np.abs(result.pressure_drop_pa - exact_dp)) / exact_dp
    )

    superficial = np.repeat(
        (result.flow_rate_m3_s / p.area_m2)[:, None], p.n_cells, axis=1
    )
    continuity_error = float(
        np.max(np.ptp(superficial, axis=1))
        / max(np.max(np.abs(superficial)), 1e-30)
    )
    pore_velocity = superficial / result.porosity

    checks = [
        Check(
            "Uniform-bed Darcy pressure",
            darcy_error < 1e-12,
            darcy_error,
            1e-12,
            "relative error",
            "Numerical pressure drop versus the exact Darcy solution.",
        ),
        Check(
            "Superficial-velocity continuity",
            continuity_error < 1e-14,
            continuity_error,
            1e-14,
            "relative range",
            "Constant-area incompressible 1D flow must have one velocity in every cell.",
        ),
        Check(
            "Finite nonnegative pore velocity",
            bool(np.all(np.isfinite(pore_velocity)) and np.all(pore_velocity >= 0)),
            float(np.min(pore_velocity)),
            0.0,
            "m/s minimum",
        ),
    ]
    return checks, {
        "z_mm": result.z_m * 1e3,
        "superficial": superficial[-1],
        "pore": pore_velocity[-1],
    }


def darcy_forchheimer_exact_check() -> Check:
    """Compare the implemented nonlinear pressure law with its exact value."""

    p = inert_case(
        hydraulic_model="darcy_forchheimer",
        flow_rate_m3_s=5.0e-6,
        permeability_initial_m2=8.0e-13,
        forchheimer_coefficient_m_inv=1.2e5,
        basket_sieve_resistance_pa_s_m3=2.0e9,
        duration_s=0.2,
    )
    result = simulate(p)
    q = p.flow_rate_m3_s / p.area_m2
    exact = p.length_m * (
        p.viscosity_pa_s * q / p.permeability_initial_m2
        + p.fluid_density_kg_m3 * p.forchheimer_coefficient_m_inv * q**2
    ) + p.flow_rate_m3_s * p.basket_sieve_resistance_pa_s_m3
    error = relative_error(float(result.pressure_drop_pa[-1]), exact)
    return Check(
        "Uniform-bed Darcy-Forchheimer pressure",
        error < 1e-12,
        error,
        1e-12,
        "relative error",
        "Numerical pressure drop versus the exact linear-plus-quadratic law.",
    )


def darcy_forchheimer_pressure_inversion_check() -> Check:
    """Check the positive quadratic root used for prescribed pressure."""

    p = inert_case(
        control_mode="pressure",
        hydraulic_model="darcy_forchheimer",
        pressure_drop_pa=7.0e5,
        pressure_pulse_duration_s=0.0,
        permeability_initial_m2=8.0e-13,
        forchheimer_coefficient_m_inv=1.2e5,
        basket_sieve_resistance_pa_s_m3=2.0e9,
        duration_s=0.2,
    )
    result = simulate(p)
    error = relative_error(float(result.pressure_drop_pa[-1]), p.pressure_drop_pa)
    return Check(
        "Darcy-Forchheimer pressure inversion",
        error < 1e-12,
        error,
        1e-12,
        "relative error",
        "The positive scalar root must reproduce the prescribed total pressure.",
    )


def washout_escape_check() -> Check:
    """Verify that a larger d_f/d_h ratio suppresses successful remobilization."""

    common = dict(
        control_mode="flow",
        hydraulic_model="darcy",
        flow_rate_m3_s=1.0e-6,
        duration_s=1.0,
        output_interval_s=0.1,
        available_fines_initial_kg_m3=0.0,
        deposited_fines_initial_kg_m3=8.0,
        release_rate_s=0.0,
        deposition_coefficient_m_inv=0.0,
        detachment_rate_s=1.0,
        detachment_exponent=1.0,
        critical_pressure_gradient_pa_m=1.0e4,
        washout_escape_model="logistic",
        escape_ratio_midpoint=0.35,
        escape_ratio_steepness=25.0,
    )
    small = simulate(replace(ModelParameters(), fine_particle_diameter_m=10e-6, **common))
    large = simulate(replace(ModelParameters(), fine_particle_diameter_m=80e-6, **common))
    deposited_difference = float(
        np.mean(large.deposited_kg_m3[-1])
        - np.mean(small.deposited_kg_m3[-1])
    )
    return Check(
        "Size-ratio washout selectivity",
        deposited_difference > 0.0,
        deposited_difference,
        0.0,
        "kg/m3 difference",
        "Larger fines must remain more strongly trapped at equal hydraulic loading.",
    )


def conservation_and_bounds_check() -> list[Check]:
    p = ModelParameters(duration_s=20.0, output_interval_s=0.25)
    result = simulate(p)
    initial = max(result.initial_fine_mass_kg, 1e-30)
    mass_error = float(np.max(np.abs(result.mass_balance_error_kg)) / initial)
    minimum_state = float(
        min(
            np.min(result.available_kg_m3),
            np.min(result.mobile_kg_m3_bulk),
            np.min(result.deposited_kg_m3),
        )
    )
    k_upper_error = float(
        max(0.0, np.max(result.permeability_m2 / p.permeability_initial_m2) - 1.0)
    )
    return [
        Check("Global fines conservation", mass_error < 1e-10, mass_error, 1e-10, "relative error"),
        Check("Nonnegative fines states", minimum_state >= -1e-12, minimum_state, -1e-12, "kg/m3 minimum"),
        Check("Permeability does not exceed K0", k_upper_error < 1e-12, k_upper_error, 1e-12, "excess ratio"),
        Check("Porosity lower bound", float(np.min(result.porosity)) >= p.minimum_porosity, float(np.min(result.porosity)), p.minimum_porosity, "minimum"),
    ]


def analytical_release_check() -> tuple[Check, dict[str, np.ndarray]]:
    p = inert_case(
        control_mode="pressure",
        pressure_drop_pa=0.0,
        duration_s=5.0,
        output_interval_s=0.1,
        available_fines_initial_kg_m3=12.0,
        release_rate_s=0.17,
        release_flow_exponent=0.0,
    )
    result = simulate(p)
    numerical = np.mean(result.available_kg_m3, axis=1)
    exact = p.available_fines_initial_kg_m3 * np.exp(-p.release_rate_s * result.time_s)
    error = float(np.max(np.abs(numerical - exact)) / p.available_fines_initial_kg_m3)
    return (
        Check("Analytical release ODE", error < 3e-3, error, 3e-3, "maximum normalized error"),
        {"time": result.time_s, "numerical": numerical, "exact": exact},
    )


def pure_advection_tvd_check() -> Check:
    p = inert_case(
        n_cells=80,
        duration_s=3.0,
        output_interval_s=0.05,
        mobile_concentration_initial_kg_m3_liquid=1.0,
    )
    result = simulate(p)
    concentration = result.concentration_kg_m3_liquid
    # Include the clean-water inlet and zero exterior outlet states. Without
    # these boundary jumps, a uniform initial field has TV=0 and the physical
    # clean-water front would be mistaken for a newly created oscillation.
    extended = np.pad(concentration, ((0, 0), (1, 1)), mode="constant")
    tv = np.sum(np.abs(np.diff(extended, axis=1)), axis=1)
    increase = float(max(0.0, np.max(tv - tv[0])))
    minimum = float(np.min(concentration))
    passed = increase < 1e-10 and minimum >= -1e-12
    return Check(
        "Pure-advection boundedness/TVD",
        passed,
        max(increase, -minimum),
        1e-10,
        "absolute",
        "First-order upwind transport should create no new oscillatory extrema.",
    )


def convergence_case(n_cells: int, max_dt: float) -> ModelParameters:
    return ModelParameters(
        control_mode="pressure",
        pressure_drop_pa=9e5,
        pressure_pulse_duration_s=0.0,
        n_cells=n_cells,
        duration_s=8.0,
        output_interval_s=0.25,
        max_time_step_s=max_dt,
        detachment_rate_s=0.0,
    )


def metrics(result) -> dict[str, float]:
    p = result.parameters
    return {
        "final_flow_ml_s": float(result.flow_rate_m3_s[-1] * 1e6),
        "outlet_fines_g": float(result.cumulative_outlet_mass_kg[-1] * 1e3),
        "minimum_K_over_K0": float(np.min(result.permeability_m2[-1] / p.permeability_initial_m2)),
        "maximum_deposit_kg_m3": float(np.max(result.deposited_kg_m3[-1])),
    }


def refinement_study(kind: str) -> pd.DataFrame:
    rows: list[dict[str, float]] = []
    if kind == "space":
        settings = [(20, 0.01), (40, 0.01), (80, 0.01), (160, 0.01)]
    else:
        settings = [(80, dt) for dt in (0.08, 0.04, 0.02, 0.01)]

    for n_cells, dt in settings:
        row = {"n_cells": n_cells, "max_time_step_s": dt}
        row.update(metrics(simulate(convergence_case(n_cells, dt))))
        rows.append(row)

    frame = pd.DataFrame(rows)
    for column in list(metrics(simulate(convergence_case(20, 0.01))).keys()):
        reference = float(frame[column].iloc[-1])
        frame[f"{column}_relative_to_finest"] = (
            np.abs(frame[column] - reference) / max(abs(reference), 1e-15)
        )
    return frame


def convergence_check(frame: pd.DataFrame, kind: str) -> Check:
    error_columns = [c for c in frame.columns if c.endswith("_relative_to_finest")]
    penultimate = float(frame[error_columns].iloc[-2].max())
    coarse = float(frame[error_columns].iloc[0].max())
    passed = penultimate < 0.08 and penultimate <= coarse + 1e-12
    return Check(
        f"{kind.capitalize()} refinement convergence",
        passed,
        penultimate,
        0.08,
        "max relative difference",
        "Difference of the penultimate resolution from the finest computed result.",
    )


def controller_checks() -> list[Check]:
    resistive = simulate(
        inert_case(control_mode="pi", duration_s=15.0, output_interval_s=0.1)
    )
    target_error_bar = abs(
        resistive.pressure_drop_pa[-1] - resistive.parameters.pressure_setpoint_pa
    ) / 1e5

    coarse = simulate(
        inert_case(
            control_mode="pi",
            permeability_initial_m2=1e-12,
            duration_s=12.0,
            output_interval_s=0.1,
        )
    )
    cap_excess = float(
        max(0.0, np.max(coarse.flow_rate_m3_s) - coarse.parameters.pump_flow_max_m3_s)
        * 1e6
    )
    return [
        Check("PI reaches attainable pressure", target_error_bar < 0.1, target_error_bar, 0.1, "bar error"),
        Check("PI respects pump-flow cap", cap_excess < 1e-10, cap_excess, 1e-10, "mL/s excess"),
        Check("Coarse case remains below target at cap", coarse.pressure_drop_pa[-1] < coarse.parameters.pressure_setpoint_pa, float(coarse.pressure_drop_pa[-1] / 1e5), float(coarse.parameters.pressure_setpoint_pa / 1e5), "bar"),
    ]


def save_plots(output_dir: Path, velocity, release, space, time) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), constrained_layout=True)
    axes[0, 0].plot(velocity["z_mm"], velocity["superficial"] * 1e3, label="Superficial")
    axes[0, 0].plot(velocity["z_mm"], velocity["pore"] * 1e3, "--", label="Pore")
    axes[0, 0].set(xlabel="Depth [mm]", ylabel="Velocity [mm/s]", title="Velocity continuity benchmark")
    axes[0, 0].legend()

    axes[0, 1].plot(release["time"], release["exact"], label="Analytical", lw=2)
    axes[0, 1].plot(release["time"], release["numerical"], "--", label="Numerical")
    axes[0, 1].set(xlabel="Time [s]", ylabel="Available fines [kg/m3]", title="Release-only benchmark")
    axes[0, 1].legend()

    error_col = "final_flow_ml_s_relative_to_finest"
    axes[1, 0].loglog(space["n_cells"].iloc[:-1], space[error_col].iloc[:-1], "o-")
    axes[1, 0].set(xlabel="Number of cells", ylabel="Relative difference", title="Spatial refinement: final flow")

    axes[1, 1].loglog(time["max_time_step_s"].iloc[:-1], time[error_col].iloc[:-1], "o-")
    axes[1, 1].invert_xaxis()
    axes[1, 1].set(xlabel="Maximum timestep [s]", ylabel="Relative difference", title="Temporal refinement: final flow")
    fig.suptitle("FineFlow-Espresso numerical verification")
    fig.savefig(output_dir / "verification_plots.png", dpi=180)
    plt.close(fig)


def write_report(output_dir: Path, checks: list[Check]) -> None:
    passed = sum(check.passed for check in checks)
    lines = [
        "FINEFLOW-ESPRESSO NUMERICAL VERIFICATION REPORT",
        "=" * 47,
        "",
        f"Overall: {passed}/{len(checks)} checks passed",
        "",
        "This report verifies selected numerical properties of the implementation.",
        "It does not validate the physical model or parameter values against data.",
        "",
    ]
    for check in checks:
        status = "PASS" if check.passed else "FAIL"
        lines.append(f"[{status}] {check.name}")
        lines.append(f"       value = {check.value:.8g} {check.units}")
        lines.append(f"       criterion/reference = {check.limit:.8g} {check.units}")
        if check.explanation:
            lines.append(f"       {check.explanation}")
        lines.append("")
    (output_dir / "verification_report.txt").write_text("\n".join(lines), encoding="utf-8")
    pd.DataFrame([check.__dict__ for check in checks]).to_csv(
        output_dir / "verification_checks.csv", index=False
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("verification_results"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    checks, velocity = hydraulic_and_velocity_checks()
    checks.append(darcy_forchheimer_exact_check())
    checks.append(darcy_forchheimer_pressure_inversion_check())
    checks.append(washout_escape_check())
    checks.extend(conservation_and_bounds_check())
    release_check, release = analytical_release_check()
    checks.append(release_check)
    checks.append(pure_advection_tvd_check())

    space = refinement_study("space")
    time = refinement_study("time")
    space.to_csv(args.output_dir / "spatial_convergence.csv", index=False)
    time.to_csv(args.output_dir / "temporal_convergence.csv", index=False)
    checks.append(convergence_check(space, "spatial"))
    checks.append(convergence_check(time, "temporal"))
    checks.extend(controller_checks())

    save_plots(args.output_dir, velocity, release, space, time)
    write_report(args.output_dir, checks)
    metadata = {
        "checks_passed": int(sum(bool(check.passed) for check in checks)),
        "checks_total": len(checks),
        "all_passed": bool(all(bool(check.passed) for check in checks)),
    }
    (args.output_dir / "verification_summary.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2))
    print(f"Verification outputs: {args.output_dir.resolve()}")
    raise SystemExit(0 if metadata["all_passed"] else 1)


if __name__ == "__main__":
    main()
