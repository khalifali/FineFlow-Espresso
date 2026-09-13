# FineFlow-Espresso

FineFlow-Espresso is a transparent one-dimensional demonstrator for fines
release, migration, deposition, blockage, detachment, and washout in an
espresso puck. It requires neither CFD nor a pore-network model.

The code is a starting framework. The included parameter values are chosen to
produce interpretable behavior and are **not calibrated espresso parameters**.
The default case holds a constant 9 bar across the puck plus basket.
Optional PI control targets 9 bar with a finite pump-flow capacity.

Complete model and numerical documentation is provided in
`documentation/fineflow_espresso_technical_report.tex` and its compiled PDF.
The report includes governing equations, assumptions, boundary conditions,
finite-volume discretization, PI control, basket accounting, CT/PNM-ready
closures, generated results, and numerical verification.

## Model states

Each axial finite-volume cell tracks:

- available fines `m_a` still attached to the coffee matrix [kg/m3 bulk],
- mobile fines inventory `epsilon*c` [kg/m3 bulk],
- deposited fines `m_d` [kg/m3 bulk],
- porosity `epsilon` [-], permeability `K` [m2], and pressure `p` [Pa].

The code conserves the sum of available, mobile, deposited, and discharged
fine mass. Clean water enters at `z=0`; fines leave advectively at `z=L`.

A level-1 basket model sits downstream of the puck. It adds a constant
hydraulic resistance and retains a constant fraction of the fines leaving the
puck. Outputs distinguish fines leaving the puck, retained by the sieve, and
reaching the cup. The controller acts on the combined puck-plus-sieve pressure.

The closures implemented in the concept are:

```text
release:     R_rel = k_rel * m_a * f_q(q)
deposition:  R_dep = k_dep * chi(z) * |q| * c * max(1 - m_d/m_cap, 0)
detachment:  R_det,* = k_det * m_d * max(G/G_crit - 1, 0)^b
washout:     R_wash = P_escape(d_f/d_h) * R_det,*
blocking:    K = K0 * exp(-beta * m_d/m_cap)
hydraulics:  |dp/dz| = mu*q/K + rho*beta_F*|q|*q
```

The deposit-to-permeability closure can be selected with
`permeability_model`: `exponential`, `kozeny_carman`, or `pnm_table`. The last
option linearly interpolates a deposit-versus-relative-permeability table
exported from pore-network calculations. Optional axial arrays
`initial_porosity_profile` and `initial_permeability_profile_m2` accept one
value per finite-volume cell, allowing depth-resolved CT/PNM inputs.

`hydraulic_model` selects `darcy` (default) or `darcy_forchheimer`. The latter is optional. It preserves CT/PNM-derived permeability in the viscous term while
adding an independently calibratable nonlinear inertial term. This retains more
structural information than calculating the entire resistance from the Ergun
porosity-and-effective-diameter correlation.

`G = |dp/dz|` is the local pressure gradient. Four operating modes are
available:

- `constant` (default): fixed `pressure_drop_pa`, with no pulse or pump cap;
- `pi`: the PI controller adjusts flow to target 9 bar;
- `pressure`: an ideal prescribed pressure, optionally including a pulse;
- `flow`: an ideal prescribed flow, with pressure predicted from resistance.

In PI mode, the current puck resistance supplies a feedforward flow estimate;
the PI correction removes residual pressure error. Pump lag prevents an
instantaneous flow change, and conditional integration prevents controller
windup. The command is clipped at `pump_flow_max_m3_s`. Consequently, a very
permeable coarse puck can remain below 9 bar even while the pump runs at its
maximum flow.

The illustrative default cap is 650 mL/min = 10.83 mL/s, corresponding to the
published maximum free-flow rating of a CEME/ULKA E5 vibration pump. This is a
reasonable demonstrator value, not a universal machine specification. A
measured machine-specific pump curve should ultimately replace the constant
cap: https://www.cemegroup.com/solenoid-pump/e5-60

## Run a named case

From this folder:

```bash
python3 -m pip install -r requirements.txt
python3 fineflow_espresso.py --config case_config.json --output-dir results
```

`case_config.json` contains named parameter sets. Change `selected_case` to
choose one; output filenames are prefixed with that name. For example,
`fine_puck` produces `fine_puck_timeseries.csv` and
`fine_puck_dashboard.png`. See `RUNNING_THE_SOLVER.txt` for a practical guide
to every parameter group and a recommended calibration order.

Optional command-line overrides:

```bash
python3 fineflow_espresso.py --pi-pressure-bar 9 --pump-cap-ml-s 10.83
python3 fineflow_espresso.py --config case_config.json --output-dir results
python3 fineflow_espresso.py --pressure-bar 9 --duration 60
python3 fineflow_espresso.py --pressure-bar 9 --duration 60 --no-pulse
python3 fineflow_espresso.py --flow-ml-s 1.0 --duration 60
python3 fineflow_espresso.py --config results/fine_puck_parameters.json
```

The output directory contains case-prefixed files:

- `*_timeseries.csv`: pressure, flow, washout, blockage position, and mass error;
- `*_profiles.npz`: complete time-dependent axial fields;
- `*_summary.json`: compact final metrics;
- `*_parameters.json`: exact parameters, reusable as a configuration;
- `*_dashboard.png`: hydraulic, washout, blockage, and permeability plots.
- `*_hydraulic_response`, `*_fines_washout`, `*_deposited_fines`, and
  `*_relative_permeability`: standalone manuscript figures in 600-dpi PNG and
  vector PDF.

The CSV and summary separately report total, puck, and sieve pressure drops,
as well as puck-outlet, basket-retained, and cup fines. In level 1, retained
fines do not change basket resistance and cannot detach again.

## Run the verification tests

```bash
python3 -m unittest -v test_fineflow.py
```

For the extended analytical and refinement verification suite:

```bash
python3 verify_fineflow.py --output-dir verification_results
```

It creates a pass/fail report, exact-versus-numerical plots, and spatial and
temporal convergence tables. This verifies implementation behavior; it does
not replace experimental validation of the physical model.

The tests check:

- conservative wash-through,
- conservative release,
- deposition-driven permeability and flow reduction,
- the detachment threshold,
- the exact Darcy-Forchheimer pressure law and nonlinear pressure inversion,
- reduced escape of fines with larger `d_f/d_h`,
- blockage erosion and increased washout during a pressure pulse,
- pressure increase under imposed flow,
- convergence to 9 bar for a sufficiently resistive puck under PI control,
- flow-cap saturation below 9 bar for a very permeable coarse puck.

## Parameters that require measurements

The most important quantities to calibrate are:

- release rate `release_rate_s`,
- deposition coefficient `deposition_coefficient_m_inv`,
- deposit capacity `deposit_capacity_kg_m3`,
- detachment rate `detachment_rate_s`,
- critical pressure gradient `critical_pressure_gradient_pa_m`,
- permeability blocking factor `permeability_blocking_factor`.

The detailed manual separates parameters into directly measured, CT-derived,
PNM-derived, hydraulic-experiment, machine, and dynamic-fines categories and
specifies a staged calibration procedure.

The machine-side quantities also require calibration: pump-flow cap, pump
time constant, and PI gains. For a specific espresso machine, record its
pressure and flow response without coffee or with several known hydraulic
resistances.

Time-resolved pressure and flow constrain total hydraulic resistance. Outlet
fine mass or turbidity is needed to separate release, deposition, and
detachment. Initial and final fines fractions or PSDs further constrain the
fine-mass inventory.

## Important limitation

Reverse flow and radial channeling are outside the scope of this work.
The axial location of blockage is an internal 1D model prediction. Global
pressure and flow alone cannot validate its physical location. Permeability recovery represents local deposit removal; it does not model
radial redistribution or channel geometry.

## Darcy selection and closure sensitivity

The default solver and active `fine_puck` case now use Darcy flow. The optional
Darcy–Forchheimer model remains available. Run `python compare_closures.py` to
reproduce the isolated inertia comparison and the separate Kozeny–Carman study.
See `closure_comparison/README.md` for the design, quantitative results and limits.

The two model selections are independent: `hydraulic_model` controls the
pressure-flow relation; `permeability_model` controls deposit-induced change in
`K/K0`. Exponential remains the active permeability choice. Its dimensionless
`permeability_blocking_factor` is an empirical blockage-strength parameter,
not a Kozeny factor. With beta = 5.5, the unclipped ratio at deposit capacity is
0.00409; this severe loss is imposed, not experimentally established.
The Kozeny-Carman option retains prescribed `K0` and changes only the relative
porosity law, assuming a fixed surface/geometry prefactor. A constant Kozeny
factor cancels in that ratio. See manual section 5 for equations and assumptions.

Finite-volume conservation follows from shared face fluxes, matched transfers
between fines inventories, and consistent basket/outlet accounting. It is
checked numerically and should not be confused with experimental validation.

The committed `results/` files now describe the constant-pressure fine-puck baseline.
The frozen `closure_comparison/` study retains its original PI settings. The technical report now
uses the controlled study in `closure_comparison/`. Manual section 12 gives
reproduction and LaTeX build commands.


## Constant pressure versus PI

```bash
# Ideal fixed total pressure: no pulse and no pump-flow cap
python3 fineflow_espresso.py --config case_config.json --pressure-bar 9 --case-name fine_constant
# Same puck, with a capped flow actuator and PI pressure target
python3 fineflow_espresso.py --config case_config.json --pi-pressure-bar 9 --case-name fine_pi
# Matched fine/coarse comparison with complete result and parameter snapshots
python3 compare_pressure_modes.py
```

`--constant-pressure-bar` is an alias for `--pressure-bar`. Both now select
`control_mode: "constant"`. Existing JSON `control_mode: "pressure"` retains
its legacy pulse behavior; use `pressure_pulse_duration_s: 0` to disable it.
Existing JSON `control_mode: "pi"` remains unchanged. A JSON file without a
mode now inherits the constant-pressure default.

The pressure is inlet minus downstream-of-basket pressure. With the outlet
reference at zero, puck-exit pressure is `Q * basket_sieve_resistance_pa_s_m3`.
Set basket resistance to zero when prescribing exactly 9 bar across the puck.
Constant mode computes Q = delta_p / (R_puck + R_basket) at each step.
It assumes the source can supply the required flow, even above the PI pump cap.
The legacy `pressure` mode supports a rectangular pulse, not an arbitrary
measured pressure trace.

Compare a fixed 9 bar run only with the corresponding 9 bar experimental
interval. For a lower, approximately steady measured pressure, use that value
with `--pressure-bar`. A changing experimental pressure requires a prescribed
history extension or a calibrated machine model. Selecting PI alone does not
make an experimental comparison valid. Pump capacity, lag and gains must be
identified independently. Saturated single-phase flow is assumed throughout;
a pressure ramp does not model initial wetting.

See `pressure_comparison/README.md` and the updated sectioned deck for the
matched comparison. The earlier PI closure comparison is retained and labeled
separately so pressure-mode and permeability-law effects are not confounded.
