"""Fit FAI/FCI at common integer durations 1--20 days to daily HBV output.

Maximum contiguous volume is searched over each complete supplied simulation.
The calendar interval displayed in the figure is not used to crop the fit data.
No event-end threshold, smoothing, normalization, or fit-quality filter is used.
"""
from pathlib import Path
from collections import defaultdict
import csv
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parent / "outputs" / "analysis"
SOURCE = ROOT / 'hydrographs.csv'
DURATIONS = np.arange(1, 21, dtype=int)
DT = 1.0
QREF, TREF = 1.0, 1.0


def write_csv(name, rows):
    with (ROOT / name).open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    before = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    groups = defaultdict(list)
    with SOURCE.open(encoding='utf-8', newline='') as stream:
        for row in csv.DictReader(stream):
            groups[row['scenario_id']].append(row)
    metrics, windows = [], []
    maximum_sum_error = 0.0
    for scenario_id, rows in groups.items():
        t = np.array([float(r['time_days']) for r in rows])
        q = np.array([float(r['Q_mm_day']) for r in rows])
        assert np.all(np.isfinite(q)) and np.all(q >= 0)
        assert np.allclose(np.diff(t), DT, rtol=0, atol=1e-12)
        assert len(q) >= int(DURATIONS[-1])
        cumulative = np.r_[0.0, np.cumsum(q * DT)]
        max_volumes, local_windows = [], []
        for duration in DURATIONS:
            # Direct summation avoids subtraction of nearly equal cumulative sums.
            sums = np.lib.stride_tricks.sliding_window_view(q, int(duration)).sum(axis=1) * DT
            independent_sums = cumulative[duration:] - cumulative[:-duration]
            err = float(np.max(np.abs(sums - independent_sums)))
            maximum_sum_error = max(maximum_sum_error, err)
            assert np.allclose(sums, independent_sums, rtol=1e-10, atol=1e-10)
            start = int(np.argmax(sums))
            volume = float(sums[start])
            assert volume > 0
            max_volumes.append(volume)
            local_windows.append(dict(panel=rows[0]['panel'], scenario_id=scenario_id,
                duration_days=int(duration), W_max_mm=volume,
                window_start_day=float(t[start] - DT/2),
                window_end_day=float(t[start + duration - 1] + DT/2)))
        x = np.log10(DURATIONS / TREF)
        y = np.log10(np.array(max_volumes) / (QREF*TREF))
        slope, intercept = np.linalg.lstsq(np.column_stack([x, np.ones(len(x))]), y, rcond=None)[0]
        predicted = intercept + slope*x
        residual = y - predicted
        r_squared = float(1 - np.sum(residual**2) / np.sum((y - y.mean())**2))
        for i, record in enumerate(local_windows):
            record.update(log10_T_over_Tref=float(x[i]), log10_W_over_Qref_Tref=float(y[i]),
                fitted_log10_W=float(predicted[i]), residual_log10_W=float(residual[i]))
        windows.extend(local_windows)
        metrics.append(dict(panel=rows[0]['panel'], scenario_id=scenario_id, label=rows[0]['label'],
            FAI=float(intercept), FCI=float(1-slope), fitted_volume_exponent=float(slope),
            R_squared=r_squared, RMSE_log10_W=float(np.sqrt(np.mean(residual**2))),
            fit_min_duration_days=1, fit_max_duration_days=20, fitted_duration_count=20,
            source_daily_cells=len(q), search_start_day=float(t[0]-DT/2),
            search_end_day=float(t[-1]+DT/2), Qref_mm_day=QREF, Tref_days=TREF,
            baseline=rows[0]['baseline']))
    by_id = {r['scenario_id']: r for r in metrics}
    for scenario_id in ['soil50', 'routing3', 'lag0']:
        for key in ['FAI', 'FCI', 'R_squared']:
            assert by_id[scenario_id][key] == by_id['baseline_a'][key]
    after = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    assert before == after
    write_csv('fitted_indices.csv', metrics)
    write_csv('maximum_window_volumes.csv', windows)
    audit = dict(source_sha256=before, source_unchanged=True, plotted_curve_count=len(metrics),
        definition='log10[W/(Qref*Tref)] = FAI + (1-FCI)*log10(T/Tref)',
        fit='Unweighted ordinary least squares; base-10 logarithms',
        fit_durations_days=DURATIONS.tolist(), Qref_mm_day=QREF, Tref_days=TREF,
        maximum_window_search='Complete daily simulation for each scenario, including routing tails',
        displayed_calendar_interval_days=[0, 20], displayed_crop_used_for_fitting=False,
        interpretation='Indices fitted at durations 1-20 days; not full-event estimates',
        event_end_threshold_used=False, fit_quality_filter_used=False, all_scenarios_retained=True,
        original_daily_values_used=True, interpolated_values_used=False,
        R_squared_min=min(r['R_squared'] for r in metrics),
        R_squared_max=max(r['R_squared'] for r in metrics),
        direct_vs_cumulative_window_sum_max_error_mm=maximum_sum_error,
        four_baseline_indices_identical=True, checks_passed=True)
    (ROOT/'index_fit_validation.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
