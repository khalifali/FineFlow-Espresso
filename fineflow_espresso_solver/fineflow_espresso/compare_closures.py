"""Controlled hydraulic-law and permeability-closure comparisons.

Run: python compare_closures.py
The original archived configuration is frozen in comparison_reference.json.
Only the hydraulic/permeability model switches differ between the three runs.
"""
from dataclasses import asdict, replace
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fineflow_espresso import _parameters_from_json, simulate

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'closure_comparison'


def compare(reference, candidate):
    """Maximum absolute and reference-peak-normalized errors, on saved times.

    Peak normalization avoids division by zero at startup. Profile errors
    include all saved times and cells, not just the final profile.
    """
    np.testing.assert_allclose(reference.time_s, candidate.time_s, rtol=0, atol=1e-12)
    metrics = {}
    for name in ['pressure_drop_pa', 'flow_rate_m3_s', 'cumulative_cup_mass_kg',
                 'deposited_kg_m3', 'permeability_m2']:
        a, b = getattr(reference, name), getattr(candidate, name)
        error = float(np.max(np.abs(a-b)))
        scale = float(np.max(np.abs(a)))
        metrics[name] = {'max_absolute_difference': error,
                        'max_difference_over_reference_peak': error / scale if scale else None}
    return metrics


def main():
    OUT.mkdir(exist_ok=True)
    _, base = _parameters_from_json(ROOT / 'comparison_reference.json')
    configurations = {
        'forchheimer_exponential': replace(base, hydraulic_model='darcy_forchheimer', permeability_model='exponential'),
        'darcy_exponential': replace(base, hydraulic_model='darcy', permeability_model='exponential'),
        'darcy_kozeny_carman': replace(base, hydraulic_model='darcy', permeability_model='kozeny_carman'),
    }
    runs = {}
    for name, parameters in configurations.items():
        print('Running', name, flush=True)
        result = simulate(parameters)
        assert result.summary()['maximum_relative_mass_balance_error'] < 1e-10
        assert np.min(result.mobile_kg_m3_bulk) >= 0
        result.save(OUT / name, case_name=name)
        runs[name] = result
    # Verify the comparison design, including all inherited dataclass defaults.
    for name, allowed in [('darcy_exponential', {'hydraulic_model'}),
                          ('darcy_kozeny_carman', {'permeability_model'})]:
        parent = 'forchheimer_exponential' if name == 'darcy_exponential' else 'darcy_exponential'
        a, b = asdict(configurations[parent]), asdict(configurations[name])
        assert {k for k in a if a[k] != b[k]} == allowed
    metrics = {
        'normalization': 'max absolute difference divided by maximum absolute reference value across all saved times/cells',
        'darcy_vs_forchheimer': compare(runs['forchheimer_exponential'], runs['darcy_exponential']),
        'kozeny_carman_vs_exponential': compare(runs['darcy_exponential'], runs['darcy_kozeny_carman']),
        'summaries': {name: r.summary() for name, r in runs.items()},
    }
    (OUT / 'comparison_metrics.json').write_text(json.dumps(metrics, indent=2)+'\n')
    pd.DataFrame(metrics['summaries']).T.to_csv(OUT / 'comparison_summary.csv', index_label='case')
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42, 'savefig.dpi': 300})
    specs = [('flow', 'flow_rate_m3_s', 1e6, 'Flow rate [mL/s]'),
             ('pressure', 'pressure_drop_pa', 1e-5, 'Total pressure drop [bar]'),
             ('cup_fines', 'cumulative_cup_mass_kg', 1e3, 'Cumulative cup fines [g]')]
    styles = [('forchheimer_exponential', 'Forchheimer + exponential', '#666666', '-'),
              ('darcy_exponential', 'Darcy + exponential', '#0065BD', '--'),
              ('darcy_kozeny_carman', 'Darcy + Kozeny–Carman', '#D55E00', '-')]
    for stem, attr, factor, ylabel in specs:
        fig, ax = plt.subplots(figsize=(6.4,4.1), layout='constrained')
        for name, label, color, ls in styles:
            r = runs[name]; ax.plot(r.time_s, getattr(r, attr)*factor, label=label, color=color, ls=ls, lw=2)
        ax.set(xlabel='Time [s]', ylabel=ylabel); ax.grid(alpha=.18); ax.legend(fontsize=9)
        for ext in ['pdf','png']:fig.savefig(OUT / f'comparison_{stem}.{ext}')
        plt.close(fig)
    for stem, attr, ylabel in [('deposit','deposited_kg_m3','Final deposited fines [kg/m³ bulk]'),
                                ('permeability','permeability_m2','Final permeability / initial permeability')]:
        fig, ax = plt.subplots(figsize=(6.4,4.1), layout='constrained')
        for name, label, color, ls in styles:
            r = runs[name]; values=getattr(r,attr)[-1]
            if stem == 'permeability': values = values/r.permeability_m2[0]
            ax.plot(r.z_m*1e3,values,label=label,color=color,ls=ls,lw=2)
        ax.set(xlabel='Axial position from inlet [mm]', ylabel=ylabel); ax.grid(alpha=.18); ax.legend(fontsize=9)
        for ext in ['pdf','png']:fig.savefig(OUT / f'comparison_{stem}.{ext}')
        plt.close(fig)
    print(json.dumps(metrics, indent=2))

if __name__ == '__main__':
    main()
