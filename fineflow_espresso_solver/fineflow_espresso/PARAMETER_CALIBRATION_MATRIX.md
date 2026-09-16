# FineFlow-Espresso numerical parameter calibration matrix

This table covers every numerical `ModelParameters` input in `fineflow_espresso.py`. Categorical switches are listed below.

| Physical name | Python syntax | Meaning / role | Intended source or fitting method |
|---|---|---|---|
| Puck height | `length_m` | Axial bed height; sets bulk volume and cell size. | **experiments-geometry** |
| Puck diameter | `diameter_m` | Sets cross-sectional area and bulk volume. | **experiments-geometry** |
| Water viscosity | `viscosity_pa_s` | Darcy viscous resistance and apparent-permeability conversion. | **experiments-fluid-temperature / property table** |
| Water density | `fluid_density_kg_m3` | Used by optional Forchheimer inertia. | **experiments-fluid-temperature / property table** |
| Number of axial cells | `n_cells` | Spatial discretization of the 1D model. | **numerical-convergence** |
| Uniform initial porosity | `porosity_initial` | Uniform epsilon0 when initial_porosity_mode=given; general baseline 0.55. | **CT-scan** |
| Coffee dose | `coffee_mass_kg` | Required for mass_height_density porosity calculation. | **experiments-dose** |
| Coffee particle/envelope density | `particle_density_kg_m3` | Required for mass_height_density porosity calculation. | **experiments-density / pycnometry** |
| Initial axial porosity profile | `initial_porosity_profile` | Optional cellwise epsilon0(z); overrides scalar in given mode. | **CT-scan** |
| Initial permeability | `permeability_initial_m2` | Clean/initial hydraulic permeability K0. | **experiments-hydraulic; CFD-DEM case01 support** |
| Initial axial permeability profile | `initial_permeability_profile_m2` | Optional K0(z) profile. | **CFD-DEM case01; experiments-hydraulic global scaling** |
| Effective deposit density | `deposit_density_kg_m3` | Converts deposited fine loading into porosity loss. | **CFD-DEM case02/case04 + CT-scan structural change** |
| Minimum porosity | `minimum_porosity` | Safety lower bound for epsilon. | **numerical/physical bound, not fitted** |
| Minimum permeability ratio | `minimum_permeability_ratio` | Regularization floor for K/K0. | **numerical-convergence / robustness** |
| Initially available fines inventory | `available_fines_initial_kg_m3` | Releasable fine mass per bulk puck volume. | **experiments-PSD + dose + experiments-massincup** |
| Initially deposited fines | `deposited_fines_initial_kg_m3` | Deposit inventory at t=0; normally zero for dry-start baseline. | **initial-condition; CFD-DEM handoff if preloaded** |
| Initial mobile-fines concentration | `mobile_concentration_initial_kg_m3_liquid` | Suspended fines in pore liquid at t=0. | **initial-condition / experiments-fluid sample** |
| Deposit capacity | `deposit_capacity_kg_m3` | Local saturation scale suppressing further capture. | **CFD-DEM case04; loaded experiments-hydraulic + mass balance** |
| Release rate constant | `release_rate_s` | Timescale for release from initially attached/releasable fines. | **experiments-massincup time history** |
| Release flow exponent | `release_flow_exponent` | Sensitivity of release rate to superficial velocity. | **experiments-massincup at several flow/pressure conditions** |
| Release reference velocity | `release_reference_velocity_m_s` | Reference q used to nondimensionalize release law. | **prescribed from experiments-hydraulic characteristic flow** |
| Axial dispersion coefficient | `axial_dispersion_m2_s` | Spreads the mobile-fines front axially. | **CFD-DEM case03/case04; transport experiment if available** |
| Deposition coefficient | `deposition_coefficient_m_inv` | Capture probability per travelled length in unsaturated pores. | **CFD-DEM case03, refined with case04 and experiments-massincup** |
| Outlet capture multiplier | `outlet_capture_multiplier` | Extra phenomenological capture weighting near outlet. | **CT-scan axial redistribution; otherwise fix to 1** |
| Outlet capture decay length | `outlet_capture_decay_length_m` | Axial length scale of outlet-weighted capture. | **CT-scan axial redistribution; otherwise not fitted** |
| Detachment rate constant | `detachment_rate_s` | Rate scale for hydraulically activated remobilization. | **CFD-DEM case05 + experiments-massincup under pressure/flow changes** |
| Critical pressure gradient | `critical_pressure_gradient_pa_m` | Threshold Gcrit below which detachment is zero. | **CFD-DEM case05 + experiments-hydraulic/massincup pressure-step tests** |
| Detachment exponent | `detachment_exponent` | Nonlinearity above Gcrit. | **CFD-DEM case05 + pressure-step experiments** |
| Exponential blocking factor | `permeability_blocking_factor` | Controls K/K0 loss versus normalized deposited loading. | **CFD-DEM case02/case04 + loaded experiments-hydraulic** |
| Deposit grid for tabulated blockage | `pnm_deposit_kg_m3` | Legacy name: deposited-loading grid for table closure. | **CFD-DEM case02/case04; loaded experiments-hydraulic** |
| Relative-permeability table | `pnm_permeability_ratio` | Legacy name: K/K0 versus deposited loading. | **CFD-DEM case02/case04 + experiments-hydraulic** |
| Forchheimer coefficient | `forchheimer_coefficient_m_inv` | Quadratic inertial pressure-loss coefficient beta_F. | **experiments-hydraulic at multiple Q; CFD-DEM case01 support** |
| Initial axial Forchheimer profile | `initial_forchheimer_profile_m_inv` | Optional beta_F,0(z). | **CFD-DEM case01; experiments-hydraulic global fit** |
| Forchheimer-permeability exponent | `forchheimer_permeability_exponent` | Optional beta_F evolution with K. | **CFD-DEM case02/case04 or loaded experiments-hydraulic; default 0** |
| Representative fine diameter | `fine_particle_diameter_m` | Effective transported fine size in accessibility/escape. | **experiments-PSD** |
| Representative coarse diameter | `coarse_effective_diameter_m` | Effective load-bearing size for hydraulic-diameter estimate. | **experiments-PSD** |
| Particle sphericity | `particle_sphericity` | Shape factor in hydraulic-diameter estimate. | **experiments-imaging / microscopy** |
| Initial hydraulic-diameter profile | `initial_hydraulic_diameter_profile_m` | Optional d_h,0(z) for geometric escape. | **CT-scan; CFD-DEM geometry check** |
| Hydraulic-diameter ratio table | `pnm_hydraulic_diameter_ratio` | Legacy name: d_h/d_h0 versus deposit loading. | **CFD-DEM case02/case04 + CT-scan structural change** |
| Escape midpoint | `escape_ratio_midpoint` | d_f/d_h ratio at 50% logistic escape probability. | **CFD-DEM case05 / fully resolved particle-passage tests** |
| Escape steepness | `escape_ratio_steepness` | Sharpness of logistic escape transition. | **CFD-DEM case05 / fully resolved particle-passage tests** |
| Size-ratio grid for tabulated escape | `pnm_size_ratio` | Legacy name: d_f/d_h grid for tabulated escape. | **CFD-DEM case05 / particle-passage campaign** |
| Escape-probability table | `pnm_escape_probability` | Legacy name: P_escape versus d_f/d_h. | **CFD-DEM case05 / particle-passage campaign** |
| Clean basket hydraulic resistance | `basket_sieve_resistance_pa_s_m3` | Pressure loss of basket/sieve external to puck. | **experiments-hydraulic empty-basket test** |
| Basket fines retention fraction | `basket_fines_retention_fraction` | Fraction of puck-outlet fines retained before cup. | **experiments-massincup + separately measured basket-retained mass** |
| Imposed pressure drop | `pressure_drop_pa` | Operating pressure for constant/pressure modes. | **experiments-hydraulic operating history; prescribed** |
| Pressure-pulse start time | `pressure_pulse_start_s` | Start of legacy ideal pressure pulse. | **experiments-operating protocol; prescribed** |
| Pressure-pulse duration | `pressure_pulse_duration_s` | Duration of legacy pressure pulse. | **experiments-operating protocol; prescribed** |
| Pressure-pulse magnitude | `pressure_pulse_pa` | Pressure during legacy pulse. | **experiments-hydraulic operating protocol; prescribed** |
| Imposed flow rate | `flow_rate_m3_s` | Operating Q for ideal flow-control mode. | **experiments-hydraulic operating condition; prescribed** |
| PI pressure setpoint | `pressure_setpoint_pa` | Machine controller target pressure. | **experiments-machine** |
| PI setpoint ramp time | `pressure_setpoint_ramp_s` | Startup ramp applied to target pressure. | **experiments-machine pressure trace** |
| Pump maximum flow | `pump_flow_max_m3_s` | Upper flow capacity of machine model. | **experiments-machine pump curve** |
| Pump time constant | `pump_time_constant_s` | First-order actuator lag. | **experiments-machine transient response** |
| PI proportional gain | `pi_kp_ml_s_per_bar` | Proportional controller gain. | **experiments-machine/controller identification** |
| PI integral gain | `pi_ki_ml_s_per_bar_s` | Integral controller gain. | **experiments-machine/controller identification** |
| Initial controller flow | `controller_initial_flow_m3_s` | Initial actuator state in PI mode. | **experiments-machine initial condition** |
| Simulation duration | `duration_s` | End time of numerical run. | **experiments-protocol / prescribed** |
| Output interval | `output_interval_s` | Saved-result sampling interval. | **numerical/output choice** |
| Maximum time step | `max_time_step_s` | Upper bound for integration step. | **numerical-convergence** |
| CFL safety factor | `cfl_safety` | Safety factor in explicit transport timestep restriction. | **numerical-convergence** |

## Categorical switches

- `initial_porosity_mode`: `given` uses CT-scan/experimental porosity; `mass_height_density` calculates epsilon0 from coffee mass, particle density, puck diameter and height. `coffee_mass_kg` and `particle_density_kg_m3` are mandatory in the latter mode.
- `permeability_model`: `exponential`, `kozeny_carman`, or the tabulated legacy `pnm_table` interface. The `pnm_` prefix is historical; tables can be populated from CFD-DEM and loaded hydraulic evidence.
- `hydraulic_model`: `darcy` baseline; `darcy_forchheimer` only when multi-flow hydraulic data justify inertia.
- `washout_escape_model`: `unity`, `logistic`, or tabulated legacy `pnm_table`.
- `control_mode`: `constant`, `pi`, `pressure`, or `flow`; operating-condition choices, not material parameters.

## CFD-DEM case mapping

- **Case 01 clean reference:** K0 and optional inertial-resistance support.
- **Case 02 held fines blockage:** K/K0 versus known retained loading; blocking law/table.
- **Case 03 dilute capture:** deposition coefficient and axial dispersion.
- **Case 04 progressive loading:** coupled retained loading, permeability and outlet mass; capacity/blockage evolution.
- **Case 05 remobilization:** detachment threshold/rate/exponent and escape under changed hydraulic forcing.

Primary calibration should still use the real coffee-puck hydraulic and fines measurements; CFD-DEM supplies mechanism-resolved support rather than replacing experimental validation.
