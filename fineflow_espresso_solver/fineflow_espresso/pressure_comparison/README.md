# Constant pressure and PI comparison

Run `python3 compare_pressure_modes.py` from the solver directory.
Each puck pair differs only in `control_mode`; exact parameters are saved.
Both use Darcy and exponential permeability. Fine puck: 60 s; coarse: 20 s.
All parameters are illustrative and uncalibrated.

| Case | Final pressure [bar] | Final flow [mL/s] | Cup fines [g] |
|---|---:|---:|---:|
| fine_puck_constant | 9.00000 | 0.044705 | 0.040234 |
| fine_puck_pi | 8.99061 | 0.045355 | 0.039930 |
| coarse_puck_constant | 9.00000 | 7.431944 | 0.051474 |
| coarse_puck_pi | 9.02217 | 7.776350 | 0.049930 |

Constant mode imposes 9 bar across puck + basket at all times. PI has a 2 s
setpoint ramp, pump lag, conditional integration and a 10.8333 mL/s cap.
The initial coarse puck requires 105.020 mL/s under ideal 9 bar. PI gives
1.658 bar at 2 s while saturated. Deposition increases resistance, and this
particular coarse case eventually reaches its target. It must not be described
as remaining below 9 bar for the whole extraction. A separate existing unit
test verifies persistent saturation for a coarse puck without evolving fines.

Maximum relative fines-mass error in these four runs is below 2.1e-15.
The tests also check the exact Darcy resistance for a nonuniform puck and
constant pressure during deposition. Numerical verification is not physical
validation. The large ideal coarse-puck flow is a consequence of the boundary
assumption, not a machine prediction or validation of Darcy at that flow.

`results/` now contains the constant-pressure fine-puck baseline.
`closure_comparison/` remains the frozen PI study. Do not compare those results
against the constant baseline as if permeability were the only changed input.
The sectioned deck includes the updated baseline and both pressure comparisons.
