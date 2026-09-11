# Controlled closure comparison

Run `python compare_closures.py` from the solver directory. The script reads the
frozen `comparison_reference.json`, runs all three cases, verifies that only the
intended switches differ, and saves parameter snapshots, time histories, profiles,
figures and comparison metrics. The active fine-puck case now uses Darcy;
Forchheimer remains selectable. The coarse-puck configuration is unchanged.

1. Forchheimer + exponential: original illustrative case.
2. Darcy + exponential: isolates removal of inertia.
3. Darcy + Kozeny–Carman: isolates the permeability closure relative to case 2.

All cases use the same initial permeability, kinetic parameters, PI controller,
40 cells, 60 s duration and timestep controls. No parameters were fitted.
Kozeny–Carman changes the relative permeability law while retaining the same K0;
this does not independently predict absolute K0 from a measured surface area.

The maximum Darcy/Forchheimer flow difference is 2.9642e-7 mL/s
(2.9102e-7 of the reference peak). Maximum pressure difference is 0.02034 Pa.
Darcy is selected for this tested regime. At 60 s, Darcy/exponential gives
0.04535 mL/s and 0.03993 g cup fines. Darcy/Kozeny–Carman gives 0.89243 mL/s
and 0.34664 g. Minimum K/K0 is 0.01263 versus 0.72264.

The permeability closure strongly affects these illustrative predictions.
This is sensitivity evidence, not experimental validation or proof that either
closure is physically correct. The original exponential law's blockage must not
be used as a target to tune Kozeny–Carman. Differences are evaluated on matching
saved times (all cells for profiles); normalization uses the reference peak, not
pointwise relative error near zero. Detailed definitions are in comparison_metrics.json.

The single-axis comparison figures are saved as vector PDF and 300 dpi PNG.
