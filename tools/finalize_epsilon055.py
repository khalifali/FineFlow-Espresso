#!/usr/bin/env python3
from pathlib import Path
import json, re, shutil
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SDIR = ROOT / 'fineflow_espresso_solver/fineflow_espresso'
DECK = ROOT / 'fineflow_espresso_sectioned_source'


def set_baseline_and_docs():
    solver = SDIR / 'fineflow_espresso.py'
    text = solver.read_text(encoding='utf-8')
    text = text.replace('porosity_initial: float = 0.36', 'porosity_initial: float = 0.55', 1)
    solver.write_text(text, encoding='utf-8')

    cfg_path = SDIR / 'case_config.json'
    cfg = json.loads(cfg_path.read_text())
    fine = cfg['cases']['fine_puck']
    fine.update(initial_porosity_mode='given', porosity_initial=0.55,
                coffee_mass_kg=None, particle_density_kg_m3=None)
    cfg_path.write_text(json.dumps(cfg, indent=2) + '\n')

    ref_path = SDIR / 'comparison_reference.json'
    ref = json.loads(ref_path.read_text())
    p = ref['parameters']
    p.update(initial_porosity_mode='given', porosity_initial=0.55,
             coffee_mass_kg=None, particle_density_kg_m3=None)
    ref_path.write_text(json.dumps(ref, indent=2) + '\n')

    readme_path = SDIR / 'README.md'
    readme = readme_path.read_text(encoding='utf-8')
    readme = readme.replace(
        '- `given`: use `porosity_initial` directly. The active `fine_puck` case currently\n  uses `porosity_initial = 0.55` as an experimentally supplied starting value.',
        '- `given`: use `porosity_initial` directly. The general solver baseline and the active\n  `fine_puck` case use `porosity_initial = 0.55` as the current CT-matched starting value.'
    )
    readme = readme.replace(
        '- `mass_height_density`: calculate\n  `epsilon0 = 1 - coffee_mass_kg/(particle_density_kg_m3*A*length_m)`.',
        '- `mass_height_density`: calculate\n  `epsilon0 = 1 - coffee_mass_kg/(particle_density_kg_m3*A*length_m)`. In this mode\n  both `coffee_mass_kg` and `particle_density_kg_m3` are mandatory positive inputs.'
    )
    if 'PARAMETER_CALIBRATION_MATRIX.md' not in readme:
        readme += '\n## Complete parameter calibration map\n\nSee `PARAMETER_CALIBRATION_MATRIX.md` for every numerical `ModelParameters` input, its physical meaning, and the intended experimental or CFD-DEM source.\n'
    readme_path.write_text(readme, encoding='utf-8')

    manual_path = SDIR / 'RUNNING_THE_SOLVER.txt'
    manual = manual_path.read_text(encoding='utf-8')
    manual = manual.replace(
        '    initial_porosity_mode = "given"\n        Use porosity_initial directly. The active fine_puck starts from 0.55.',
        '    initial_porosity_mode = "given"\n        Use porosity_initial directly. The general baseline and active fine_puck start from 0.55.'
    )
    manual = manual.replace(
        '    initial_porosity_mode = "mass_height_density"\n        Compute epsilon0 = 1 - m_coffee/(rho_particle*A*H), using\n        coffee_mass_kg, particle_density_kg_m3, puck area A, and H=length_m.',
        '    initial_porosity_mode = "mass_height_density"\n        Compute epsilon0 = 1 - m_coffee/(rho_particle*A*H), using puck area A\n        and H=length_m. coffee_mass_kg and particle_density_kg_m3 are REQUIRED\n        positive inputs in this mode; the solver raises an error if either is missing.'
    )
    manual_path.write_text(manual, encoding='utf-8')


def write_matrix():
    rows = [
    ('Puck height','length_m','Axial bed height; sets bulk volume and cell size.','experiments-geometry'),
    ('Puck diameter','diameter_m','Sets cross-sectional area and bulk volume.','experiments-geometry'),
    ('Water viscosity','viscosity_pa_s','Darcy viscous resistance and apparent-permeability conversion.','experiments-fluid-temperature / property table'),
    ('Water density','fluid_density_kg_m3','Used by optional Forchheimer inertia.','experiments-fluid-temperature / property table'),
    ('Number of axial cells','n_cells','Spatial discretization of the 1D model.','numerical-convergence'),
    ('Uniform initial porosity','porosity_initial','Uniform epsilon0 when initial_porosity_mode=given; general baseline 0.55.','CT-scan'),
    ('Coffee dose','coffee_mass_kg','Required for mass_height_density porosity calculation.','experiments-dose'),
    ('Coffee particle/envelope density','particle_density_kg_m3','Required for mass_height_density porosity calculation.','experiments-density / pycnometry'),
    ('Initial axial porosity profile','initial_porosity_profile','Optional cellwise epsilon0(z); overrides scalar in given mode.','CT-scan'),
    ('Initial permeability','permeability_initial_m2','Clean/initial hydraulic permeability K0.','experiments-hydraulic; CFD-DEM case01 support'),
    ('Initial axial permeability profile','initial_permeability_profile_m2','Optional K0(z) profile.','CFD-DEM case01; experiments-hydraulic global scaling'),
    ('Effective deposit density','deposit_density_kg_m3','Converts deposited fine loading into porosity loss.','CFD-DEM case02/case04 + CT-scan structural change'),
    ('Minimum porosity','minimum_porosity','Safety lower bound for epsilon.','numerical/physical bound, not fitted'),
    ('Minimum permeability ratio','minimum_permeability_ratio','Regularization floor for K/K0.','numerical-convergence / robustness'),
    ('Initially available fines inventory','available_fines_initial_kg_m3','Releasable fine mass per bulk puck volume.','experiments-PSD + dose + experiments-massincup'),
    ('Initially deposited fines','deposited_fines_initial_kg_m3','Deposit inventory at t=0; normally zero for dry-start baseline.','initial-condition; CFD-DEM handoff if preloaded'),
    ('Initial mobile-fines concentration','mobile_concentration_initial_kg_m3_liquid','Suspended fines in pore liquid at t=0.','initial-condition / experiments-fluid sample'),
    ('Deposit capacity','deposit_capacity_kg_m3','Local saturation scale suppressing further capture.','CFD-DEM case04; loaded experiments-hydraulic + mass balance'),
    ('Release rate constant','release_rate_s','Timescale for release from initially attached/releasable fines.','experiments-massincup time history'),
    ('Release flow exponent','release_flow_exponent','Sensitivity of release rate to superficial velocity.','experiments-massincup at several flow/pressure conditions'),
    ('Release reference velocity','release_reference_velocity_m_s','Reference q used to nondimensionalize release law.','prescribed from experiments-hydraulic characteristic flow'),
    ('Axial dispersion coefficient','axial_dispersion_m2_s','Spreads the mobile-fines front axially.','CFD-DEM case03/case04; transport experiment if available'),
    ('Deposition coefficient','deposition_coefficient_m_inv','Capture probability per travelled length in unsaturated pores.','CFD-DEM case03, refined with case04 and experiments-massincup'),
    ('Outlet capture multiplier','outlet_capture_multiplier','Extra phenomenological capture weighting near outlet.','CT-scan axial redistribution; otherwise fix to 1'),
    ('Outlet capture decay length','outlet_capture_decay_length_m','Axial length scale of outlet-weighted capture.','CT-scan axial redistribution; otherwise not fitted'),
    ('Detachment rate constant','detachment_rate_s','Rate scale for hydraulically activated remobilization.','CFD-DEM case05 + experiments-massincup under pressure/flow changes'),
    ('Critical pressure gradient','critical_pressure_gradient_pa_m','Threshold Gcrit below which detachment is zero.','CFD-DEM case05 + experiments-hydraulic/massincup pressure-step tests'),
    ('Detachment exponent','detachment_exponent','Nonlinearity above Gcrit.','CFD-DEM case05 + pressure-step experiments'),
    ('Exponential blocking factor','permeability_blocking_factor','Controls K/K0 loss versus normalized deposited loading.','CFD-DEM case02/case04 + loaded experiments-hydraulic'),
    ('Deposit grid for tabulated blockage','pnm_deposit_kg_m3','Legacy name: deposited-loading grid for table closure.','CFD-DEM case02/case04; loaded experiments-hydraulic'),
    ('Relative-permeability table','pnm_permeability_ratio','Legacy name: K/K0 versus deposited loading.','CFD-DEM case02/case04 + experiments-hydraulic'),
    ('Forchheimer coefficient','forchheimer_coefficient_m_inv','Quadratic inertial pressure-loss coefficient beta_F.','experiments-hydraulic at multiple Q; CFD-DEM case01 support'),
    ('Initial axial Forchheimer profile','initial_forchheimer_profile_m_inv','Optional beta_F,0(z).','CFD-DEM case01; experiments-hydraulic global fit'),
    ('Forchheimer-permeability exponent','forchheimer_permeability_exponent','Optional beta_F evolution with K.','CFD-DEM case02/case04 or loaded experiments-hydraulic; default 0'),
    ('Representative fine diameter','fine_particle_diameter_m','Effective transported fine size in accessibility/escape.','experiments-PSD'),
    ('Representative coarse diameter','coarse_effective_diameter_m','Effective load-bearing size for hydraulic-diameter estimate.','experiments-PSD'),
    ('Particle sphericity','particle_sphericity','Shape factor in hydraulic-diameter estimate.','experiments-imaging / microscopy'),
    ('Initial hydraulic-diameter profile','initial_hydraulic_diameter_profile_m','Optional d_h,0(z) for geometric escape.','CT-scan; CFD-DEM geometry check'),
    ('Hydraulic-diameter ratio table','pnm_hydraulic_diameter_ratio','Legacy name: d_h/d_h0 versus deposit loading.','CFD-DEM case02/case04 + CT-scan structural change'),
    ('Escape midpoint','escape_ratio_midpoint','d_f/d_h ratio at 50% logistic escape probability.','CFD-DEM case05 / fully resolved particle-passage tests'),
    ('Escape steepness','escape_ratio_steepness','Sharpness of logistic escape transition.','CFD-DEM case05 / fully resolved particle-passage tests'),
    ('Size-ratio grid for tabulated escape','pnm_size_ratio','Legacy name: d_f/d_h grid for tabulated escape.','CFD-DEM case05 / particle-passage campaign'),
    ('Escape-probability table','pnm_escape_probability','Legacy name: P_escape versus d_f/d_h.','CFD-DEM case05 / particle-passage campaign'),
    ('Clean basket hydraulic resistance','basket_sieve_resistance_pa_s_m3','Pressure loss of basket/sieve external to puck.','experiments-hydraulic empty-basket test'),
    ('Basket fines retention fraction','basket_fines_retention_fraction','Fraction of puck-outlet fines retained before cup.','experiments-massincup + separately measured basket-retained mass'),
    ('Imposed pressure drop','pressure_drop_pa','Operating pressure for constant/pressure modes.','experiments-hydraulic operating history; prescribed'),
    ('Pressure-pulse start time','pressure_pulse_start_s','Start of legacy ideal pressure pulse.','experiments-operating protocol; prescribed'),
    ('Pressure-pulse duration','pressure_pulse_duration_s','Duration of legacy pressure pulse.','experiments-operating protocol; prescribed'),
    ('Pressure-pulse magnitude','pressure_pulse_pa','Pressure during legacy pulse.','experiments-hydraulic operating protocol; prescribed'),
    ('Imposed flow rate','flow_rate_m3_s','Operating Q for ideal flow-control mode.','experiments-hydraulic operating condition; prescribed'),
    ('PI pressure setpoint','pressure_setpoint_pa','Machine controller target pressure.','experiments-machine'),
    ('PI setpoint ramp time','pressure_setpoint_ramp_s','Startup ramp applied to target pressure.','experiments-machine pressure trace'),
    ('Pump maximum flow','pump_flow_max_m3_s','Upper flow capacity of machine model.','experiments-machine pump curve'),
    ('Pump time constant','pump_time_constant_s','First-order actuator lag.','experiments-machine transient response'),
    ('PI proportional gain','pi_kp_ml_s_per_bar','Proportional controller gain.','experiments-machine/controller identification'),
    ('PI integral gain','pi_ki_ml_s_per_bar_s','Integral controller gain.','experiments-machine/controller identification'),
    ('Initial controller flow','controller_initial_flow_m3_s','Initial actuator state in PI mode.','experiments-machine initial condition'),
    ('Simulation duration','duration_s','End time of numerical run.','experiments-protocol / prescribed'),
    ('Output interval','output_interval_s','Saved-result sampling interval.','numerical/output choice'),
    ('Maximum time step','max_time_step_s','Upper bound for integration step.','numerical-convergence'),
    ('CFL safety factor','cfl_safety','Safety factor in explicit transport timestep restriction.','numerical-convergence'),
    ]
    out = ['# FineFlow-Espresso numerical parameter calibration matrix','',
           'This table covers every numerical `ModelParameters` input in `fineflow_espresso.py`. Categorical switches are listed below.','',
           '| Physical name | Python syntax | Meaning / role | Intended source or fitting method |','|---|---|---|---|']
    for a,b,c,d in rows:
        out.append(f'| {a} | `{b}` | {c} | **{d}** |')
    out += ['', '## Categorical switches', '',
            '- `initial_porosity_mode`: `given` uses CT-scan/experimental porosity; `mass_height_density` calculates epsilon0 from coffee mass, particle density, puck diameter and height. `coffee_mass_kg` and `particle_density_kg_m3` are mandatory in the latter mode.',
            '- `permeability_model`: `exponential`, `kozeny_carman`, or the tabulated legacy `pnm_table` interface. The `pnm_` prefix is historical; tables can be populated from CFD-DEM and loaded hydraulic evidence.',
            '- `hydraulic_model`: `darcy` baseline; `darcy_forchheimer` only when multi-flow hydraulic data justify inertia.',
            '- `washout_escape_model`: `unity`, `logistic`, or tabulated legacy `pnm_table`.',
            '- `control_mode`: `constant`, `pi`, `pressure`, or `flow`; operating-condition choices, not material parameters.', '',
            '## CFD-DEM case mapping', '',
            '- **Case 01 clean reference:** K0 and optional inertial-resistance support.',
            '- **Case 02 held fines blockage:** K/K0 versus known retained loading; blocking law/table.',
            '- **Case 03 dilute capture:** deposition coefficient and axial dispersion.',
            '- **Case 04 progressive loading:** coupled retained loading, permeability and outlet mass; capacity/blockage evolution.',
            '- **Case 05 remobilization:** detachment threshold/rate/exponent and escape under changed hydraulic forcing.', '',
            'Primary calibration should still use the real coffee-puck hydraulic and fines measurements; CFD-DEM supplies mechanism-resolved support rather than replacing experimental validation.']
    (SDIR / 'PARAMETER_CALIBRATION_MATRIX.md').write_text('\n'.join(out)+'\n', encoding='utf-8')


def sync_after_runs():
    results = SDIR / 'results'
    figdir = DECK / 'figures'
    for srcname,dstname in [
        ('fine_puck_hydraulic_response.pdf','hydraulic_response.pdf'),('fine_puck_hydraulic_response.png','hydraulic_response.png'),
        ('fine_puck_fines_washout.pdf','fines_washout.pdf'),('fine_puck_fines_washout.png','fines_washout.png'),
        ('fine_puck_deposited_fines.pdf','deposited_fines.pdf'),('fine_puck_deposited_fines.png','deposited_fines.png'),
        ('fine_puck_relative_permeability.pdf','relative_permeability.pdf'),('fine_puck_relative_permeability.png','relative_permeability.png'),
        ('fine_puck_porosity_profiles.pdf','porosity_profiles.pdf'),('fine_puck_porosity_profiles.png','porosity_profiles.png')]:
        shutil.copy2(results/srcname, figdir/dstname)
    pcomp=SDIR/'pressure_comparison'; ccomp=SDIR/'closure_comparison'
    for n in ['fine_puck_hydraulics.pdf','coarse_puck_hydraulics.pdf','coarse_puck_fines.pdf']:
        shutil.copy2(pcomp/n, figdir/n)
    for n in ['comparison_flow.pdf','comparison_permeability.pdf','comparison_cup_fines.pdf']:
        shutil.copy2(ccomp/n, figdir/n)

    base=json.loads((results/'fine_puck_summary.json').read_text())
    prof=np.load(results/'fine_puck_profiles.npz')
    emin,emax=float(prof['porosity'][-1].min()),float(prof['porosity'][-1].max())
    pressure={r['case']:r for r in json.loads((pcomp/'comparison_summary.json').read_text())}
    metrics=json.loads((ccomp/'comparison_metrics.json').read_text())
    summ=metrics['summaries']; dx=summ['darcy_exponential']; kx=summ['darcy_kozeny_carman']; fx=summ['forchheimer_exponential']
    maxdp=metrics['darcy_vs_forchheimer']['pressure_drop_pa']['max_absolute_difference']

    texpath=DECK/'fineflow_espresso_sectioned.tex'; tex=texpath.read_text(encoding='utf-8')
    tex=re.sub(r'Pressure \[bar\] & [0-9.]+ & [0-9.]+', f"Pressure [bar] & {fx['final_pressure_drop_bar']:.4f} & {dx['final_pressure_drop_bar']:.4f}", tex, count=1)
    tex=re.sub(r'Flow rate \[mL/s\] & [0-9.]+ & [0-9.]+', f"Flow rate [mL/s] & {fx['final_flow_rate_ml_s']:.5f} & {dx['final_flow_rate_ml_s']:.5f}", tex, count=1)
    tex=re.sub(r'Cup fines \[g\] & [0-9.]+ & [0-9.]+', f"Cup fines [g] & {fx['cumulative_cup_fines_g']:.5f} & {dx['cumulative_cup_fines_g']:.5f}", tex, count=1)
    tex=re.sub(r'largest pressure difference was only \\textbf\{[^}]+ Pa\}', f'largest pressure difference was only \\textbf{{{maxdp:.4f} Pa}}', tex, count=1)
    tex=re.sub(r'\\item flow: [0-9.]+ mL/s;', f"\\item flow: {base['final_flow_rate_ml_s']:.4f} mL/s;", tex, count=1)
    tex=re.sub(r'\\item minimum \$K/K_0\$: [0-9.]+\.', f"\\item minimum $K/K_0$: {base['minimum_local_permeability_ratio']:.4f}.", tex, count=1)
    tex=re.sub(r'\\item puck outlet: [0-9.]+ g;', f"\\item puck outlet: {base['cumulative_puck_outlet_fines_g']:.4f} g;", tex, count=1)
    tex=re.sub(r'\\item basket retained: [0-9.]+ g;', f"\\item basket retained: {base['cumulative_basket_retained_fines_g']:.4f} g;", tex, count=1)
    tex=re.sub(r'\\item cup: [0-9.]+ g\.', f"\\item cup: {base['cumulative_cup_fines_g']:.4f} g.", tex, count=1)

    start=tex.index('\\begin{frame}{Porosity profile evolves as fines redistribute}')
    end=tex.index('\\begin{frame}{Fine puck: constant pressure versus PI}', start)
    slide=rf'''\begin{{frame}}{{Porosity profile evolves as fines redistribute}}
{{\small\color{{tumblue}}Dry baseline $\varepsilon_0=0.55$ chosen from the CT dry-state level. Constant 9 bar across puck + basket.}}\par\medskip
\begin{{columns}}[T,onlytextwidth]
 \column{{.66\textwidth}}\centering
  \safeimage{{porosity_profiles.pdf}}{{\linewidth}}{{4.05cm}}
 \column{{.30\textwidth}}\small
  \hi{{Predicted structure}}
  \begin{{itemize}}
   \item Direct model observable: $\varepsilon(z,t)$.
   \item Dry baseline: $\varepsilon_0=0.55$.
   \item Deposition lowers porosity, strongest near outlet/bottom.
   \item At 60 s: $\varepsilon\approx{emin:.3f}$--${emax:.3f}$.
  \end{{itemize}}
\end{{columns}}
\vspace{{.08cm}}
\takeaway{{After matching the dry-state baseline, compare magnitude, axial trend and time evolution with interrupted CT profiles.}}
\slidesource{{Illustrative 1D result; the GL22, GL25 and GL28 dry CT profiles have similar mean porosity but different local structure.}}
\end{{frame}}

'''
    tex=tex[:start]+slide+tex[end:]
    fc=pressure['fine_puck_constant']; fp=pressure['fine_puck_pi']; cc=pressure['coarse_puck_constant']
    tex=re.sub(r'At 60 s: constant pressure gives \$Q=[0-9.]+\$ mL/s; PI gives \$Q=[0-9.]+\$ mL/s\.', f"At 60 s: constant pressure gives $Q={fc['final_flow_rate_ml_s']:.4f}$ mL/s; PI gives $Q={fp['final_flow_rate_ml_s']:.4f}$ mL/s.", tex, count=1)
    tex=re.sub(r'Constant 9 bar initially requires [0-9.]+ mL/s, above the PI cap of 10\.83 mL/s\. PI reaches only [0-9.]+ bar at 2 s\.', f"Constant 9 bar initially requires {cc['initial_flow_rate_ml_s']:.1f} mL/s, above the PI cap of {cc['pump_flow_cap_ml_s']:.2f} mL/s. PI is initially flow-limited before resistance increases.", tex, count=1)
    tex=re.sub(r'Exponential:\\\\$Q=[0-9.]+\$ mL/s', f"Exponential:\\\\$Q={dx['final_flow_rate_ml_s']:.3f}$ mL/s", tex, count=1)
    tex=re.sub(r'Kozeny--Carman:\\\\$Q=[0-9.]+\$ mL/s', f"Kozeny--Carman:\\\\$Q={kx['final_flow_rate_ml_s']:.3f}$ mL/s", tex, count=1)
    tex=re.sub(r'Exponential:\\\\$0\.[0-9]+\$', f"Exponential:\\\\${dx['minimum_local_permeability_ratio']:.4f}$", tex, count=1)
    tex=re.sub(r'Kozeny--Carman:\\\\$0\.[0-9]+\$', f"Kozeny--Carman:\\\\${kx['minimum_local_permeability_ratio']:.3f}$", tex, count=1)
    tex=re.sub(r'Exponential:\\\\$[0-9.]+\$ g', f"Exponential:\\\\${dx['cumulative_cup_fines_g']:.4f}$ g", tex, count=1)
    tex=re.sub(r'Kozeny--Carman:\\\\$[0-9.]+\$ g', f"Kozeny--Carman:\\\\${kx['cumulative_cup_fines_g']:.4f}$ g", tex, count=1)
    texpath.write_text(tex, encoding='utf-8')

    trpath=DECK/'transcript.json'; tr=json.loads(trpath.read_text())
    tr['32']=f"The controlled inertia test uses the common epsilon0=0.55 baseline. Darcy and Darcy-Forchheimer remain practically identical; the maximum pressure difference is {maxdp:.4f} pascals."
    tr['33']=f"With epsilon0=0.55, the constant-pressure Darcy-exponential case ends at {base['final_flow_rate_ml_s']:.4f} millilitres per second and minimum K over K0 {base['minimum_local_permeability_ratio']:.4f}. The result is illustrative and closure-dependent."
    tr['34']=f"At sixty seconds, {base['cumulative_puck_outlet_fines_g']:.4f} grams leave the puck, {base['cumulative_basket_retained_fines_g']:.4f} grams are retained by the basket, and {base['cumulative_cup_fines_g']:.4f} grams reach the cup."
    tr['36']=f"The dry baseline is epsilon0=0.55, selected from the approximate mean CT dry-state level. At sixty seconds the illustrative axial range is {emin:.3f} to {emax:.3f}. We now compare the magnitude, axial trend and time evolution rather than an artificial initial offset."
    trpath.write_text(json.dumps(tr, indent=2, ensure_ascii=False)+'\n')

    rb=DECK/'README_BUILD.txt'; txt=rb.read_text()
    txt=txt.replace('Baseline field figures now use constant 9 bar; the earlier closure comparison retains PI settings.', 'Baseline field figures use constant 9 bar and epsilon0=0.55; pressure and closure comparisons are regenerated with the 0.55 baseline.')
    rb.write_text(txt)


if __name__ == '__main__':
    import argparse
    ap=argparse.ArgumentParser(); ap.add_argument('stage', choices=['prepare','sync']); args=ap.parse_args()
    if args.stage=='prepare':
        set_baseline_and_docs(); write_matrix()
    else:
        sync_after_runs()
