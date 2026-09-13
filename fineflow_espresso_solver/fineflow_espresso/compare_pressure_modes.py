"""Run matched constant/PI studies: python3 compare_pressure_modes.py.

Outputs complete runs, snapshots, tables and slide figures in pressure_comparison/.
Both pucks use Darcy + exponential. Coarse fines remain active and may increase
resistance enough to escape saturation. Parameters are illustrative, not fitted.
"""
from dataclasses import asdict, replace
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fineflow_espresso import ModelParameters, simulate
ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'pressure_comparison'


def main():
    OUT.mkdir(exist_ok=True)
    cases = json.loads((ROOT / 'case_config.json').read_text())['cases']
    summaries = []
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False, 'pdf.fonttype': 42})
    for case in ('fine_puck', 'coarse_puck'):
        values = {k: v for k, v in cases[case].items() if not k.startswith('_')}
        base = replace(ModelParameters(), **values)
        base = replace(base, hydraulic_model='darcy', permeability_model='exponential',
                       pressure_drop_pa=9e5, pressure_setpoint_pa=9e5,
                       pressure_pulse_duration_s=0.0)
        runs = {}
        for mode in ('constant', 'pi'):
            name = f'{case}_{mode}'
            print('Running', name, flush=True)
            result = simulate(replace(base, control_mode=mode))
            result.save(OUT / name, case_name=name)
            assert result.summary()['maximum_relative_mass_balance_error'] < 1e-10
            if mode == 'constant':
                np.testing.assert_allclose(result.pressure_drop_pa, 9e5, rtol=1e-12)
            else:
                assert np.max(result.flow_rate_m3_s) <= base.pump_flow_max_m3_s * (1+1e-12)
            summaries.append({'case': name, **result.summary()})
            runs[mode] = result
        a, b = (asdict(r.parameters) for r in runs.values())
        assert {k for k in a if a[k] != b[k]} == {'control_mode'}
        specs = [('pressure_drop_pa', 1e-5, 'Total pressure drop [bar]'),
                 ('flow_rate_m3_s', 1e6, 'Flow rate [mL/s]'),
                 ('cumulative_cup_mass_kg', 1e3, 'Cumulative cup fines [g]')]
        for stem, selected in [('hydraulics', specs[:2]), ('fines', specs[2:])]:
            fig, axes = plt.subplots(1, len(selected), figsize=(10 if len(selected)==2 else 6.4, 3.7), squeeze=False, layout='constrained')
            for ax, (attr, factor, label) in zip(axes.flat, selected):
                for mode, color, style in [('constant', '#0065BD', '-'), ('pi', '#D55E00', '--')]:
                    r = runs[mode]
                    ax.plot(r.time_s, getattr(r, attr)*factor, color=color, ls=style, lw=2,
                            label='Constant 9 bar' if mode=='constant' else 'PI target 9 bar')
                if attr == 'flow_rate_m3_s' and case == 'coarse_puck':
                    ax.axhline(base.pump_flow_max_m3_s*1e6, color='0.4', ls=':', label='PI pump cap')
                ax.set(xlabel='Time [s]', ylabel=label)
                ax.grid(alpha=.2)
                ax.legend(fontsize=9)
            for ext in ('pdf', 'png'):
                fig.savefig(OUT / f'{case}_{stem}.{ext}', dpi=200)
            plt.close(fig)
    (OUT / 'comparison_summary.json').write_text(json.dumps(summaries, indent=2)+'\n')
    pd.DataFrame(summaries).to_csv(OUT / 'comparison_summary.csv', index=False)
    print(pd.DataFrame(summaries)[['case', 'final_pressure_drop_bar', 'final_flow_rate_ml_s']].to_string(index=False))


if __name__ == '__main__':
    main()
