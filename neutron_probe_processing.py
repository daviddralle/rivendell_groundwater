"""Calibrate the frozen neutron-probe observations and integrate fixed intervals.

Run by notebook 01. Source tables are never modified. Calibration coefficients,
profile coverage, exclusions, and integration bounds are exported for auditing.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import TheilSenRegressor
from analysis_utils import integrate_profile

REFERENCE_PROBE = 'H340702152_01'
REFERENCE_DATE = pd.Timestamp('2015-10-25')
# 0.75 m is below the reported soil profile. Upper bounds retain repeated
# vadose-zone survey coverage; deeper, intermittently sampled rock is excluded.
DOMAINS = {2: (0.75, 9.), 5: (0.75, 10.), 6: (0.75, 4.5),
           7: (0.75, 4.), 12: (0.75, 3.), 13: (0.75, 5.),
           14: (0.75, 5.), 15: (0.75, 10.), 16: (0.75, 10.)}


def assign_calibration_periods(observations, intervals):
    """Match the instrument identity AND its half-open calibration period."""
    matched = pd.merge_asof(
        observations.sort_values('Date'), intervals.sort_values('Start'),
        left_on='Date', right_on='Start', by='Serial Number', direction='backward')
    valid = matched['Start'].notna() & (matched['Date'] < matched['End'])
    return matched.loc[valid].copy(), matched.loc[~valid].copy()


def process_neutron_profiles(data_dir='Data'):
    base = Path(data_dir)
    raw = pd.read_csv(base/'raw/neutron_probe_rivendell_through_wy2025.csv', low_memory=False)
    assert 'source_row' in raw  # Original frozen-export row, retained by the subset.
    raw = raw[(raw['Site'] == 'Rivendell') & raw['Flag'].isna()].copy()
    raw['Date'] = pd.to_datetime(raw['Date'], errors='coerce')
    for col in ['Cable length (ft)', 'Neutron Count (MD)', 'Hole ID']:
        raw[col] = pd.to_numeric(raw[col], errors='coerce')
    raw = raw.dropna(subset=['Date', 'Hole ID', 'Cable length (ft)', 'Neutron Count (MD)'])
    raw['Hole ID'] = raw['Hole ID'].astype(int)
    raw = raw[raw['Hole ID'].isin([1,2,3,5,6,7,10,12,13,14,15,16])].copy()
    count_ok = raw['Neutron Count (MD)'].between(0, 9000)
    excluded = [raw.loc[~count_ok].assign(reason='Count outside 0–9000')]
    raw = raw.loc[count_ok].drop_duplicates(
        ['Date', 'Hole ID', 'Cable length (ft)', 'Serial Number', 'Neutron Count (MD)']).copy()
    intervals = pd.read_csv(base/'raw/probe_conversion_intervals.csv').rename(columns={
        'Date Lower':'Start', 'Date Upper':'End', 'New Conversion ID':'Conversion ID'})
    for col in ['Start', 'End']:
        intervals[col] = pd.to_datetime(intervals[col])
    for serial, group in intervals.groupby('Serial Number'):
        group = group.sort_values('Start')
        assert not (group['Start'].iloc[1:].to_numpy() < group['End'].iloc[:-1].to_numpy()).any(), serial
    cal, unmatched = assign_calibration_periods(
        raw, intervals[['Start', 'End', 'Serial Number', 'Conversion ID']])
    excluded.append(unmatched.assign(reason='No matching serial-number calibration period'))

    # Repeated readings do not increase the number of independent calibration dates.
    ports = cal.groupby(['Hole ID','Cable length (ft)','Conversion ID','Date'])['Neutron Count (MD)'].median()
    stats = ports.groupby(level=[0,1,2]).agg(['median','mean','std','count']).reset_index()
    stats['cv'] = stats['std'] / stats['mean']
    reference = stats[(stats['Conversion ID'] == REFERENCE_PROBE) &
                      (stats['cv'] < 0.02) & (stats['count'] >= 10)]
    valid = stats.merge(reference[['Hole ID','Cable length (ft)']], on=['Hole ID','Cable length (ft)'])
    pivot = valid.pivot(index=['Hole ID','Cable length (ft)'], columns='Conversion ID', values='median')
    known = {REFERENCE_PROBE: dict(slope=1., intercept=0., parent='', chain_depth=0, n_ports=len(reference), rmse_counts=0.)}
    remaining = sorted(set(cal['Conversion ID']) - set(known))
    while remaining:
        candidates = []
        for target in remaining:
            if target not in pivot:
                continue
            for parent in known:
                pairs = pivot[[target, parent]].dropna()
                if len(pairs) >= 3:
                    candidates.append((known[parent]['chain_depth'], -len(pairs), target, parent, pairs))
        if not candidates:
            break
        _, _, target, parent, pairs = min(candidates, key=lambda item: item[:4])
        fit = TheilSenRegressor(random_state=42).fit(pairs[[target]].values, pairs[parent].values)
        slope, intercept = float(fit.coef_[0]), float(fit.intercept_)
        assert slope > 0, target
        params = known[parent]
        known[target] = dict(slope=params['slope']*slope,
                             intercept=params['slope']*intercept+params['intercept'],
                             parent=parent, chain_depth=params['chain_depth']+1,
                             n_ports=len(pairs), rmse_counts=float(np.sqrt(np.mean((fit.predict(pairs[[target]].values)-pairs[parent])**2))))
        remaining.remove(target)
    excluded.append(cal[~cal['Conversion ID'].isin(known)].assign(reason='Insufficient invariant ports for calibration'))
    cal = cal[cal['Conversion ID'].isin(known)].copy()
    calibration = pd.DataFrame.from_dict(known, orient='index').rename_axis('Conversion ID')
    calibration.to_csv(base/'neutron_probe_calibration.csv')
    cal['Adjusted Count'] = cal['Neutron Count (MD)'] * cal['Conversion ID'].map(calibration.slope) + cal['Conversion ID'].map(calibration.intercept)

    casing = pd.read_csv(base/'raw/rivendell_borehole_stickup.csv')
    casing = casing[casing.Site == 'Rivendell'].copy()
    casing['Hole ID'] = pd.to_numeric(casing['Hole ID']).astype(int)
    casing['Casing diameter (in)'] = pd.to_numeric(casing['Casing diameter (in)'])
    cal = cal.merge(casing[['Hole ID','Wellhead height (m)','Casing diameter (in)']], on='Hole ID', validate='many_to_one')
    cal['Depth (m)'] = cal['Cable length (ft)']*0.3048 - cal['Wellhead height (m)'] - 0.3429
    # Original barrel calibrations, indexed by actual casing diameter.
    slope = cal['Casing diameter (in)'].map({2.:0.000048, 3.:0.000072})
    intercept = cal['Casing diameter (in)'].map({2.:-0.015400, 3.:-0.029900})
    assert slope.notna().all()
    cal['theta'] = cal['Adjusted Count']*slope + intercept
    failed = cal.groupby(['Hole ID','Date'])['theta'].transform('max') > 1
    excluded.append(cal.loc[failed].assign(reason='Calibrated theta > 1 in profile'))
    cal = cal.loc[~failed].copy()
    excluded.append(cal[cal.theta < 0].assign(reason='Negative calibrated water content'))
    cal = cal[cal.theta >= 0].copy()
    pd.concat(excluded, ignore_index=True).to_csv(base/'neutron_probe_excluded_observations.csv', index=False)
    cal.to_csv(base/'neutron_probe_calibrated_observations.csv', index=False)
    profiles = cal.groupby(['Hole ID','Date','Depth (m)'], as_index=False)['theta'].median()
    profiles.to_csv(base/'rivendell_theta_depth_timeseries.csv', index=False)

    # Restrict the manuscript to complete water years through WY2025.
    selected = profiles[profiles.Date.between('2014-01-01','2025-09-30')].copy()
    gw = pd.read_csv(base/'groundwater_hourly_resampled.csv', index_col=0, parse_dates=True).resample('D').median()
    records, coverage = [], []
    for well, (lower, upper) in DOMAINS.items():
        reference = profiles[(profiles['Hole ID'] == well) & (profiles.Date == REFERENCE_DATE)].sort_values('Depth (m)')
        ref_integral = integrate_profile(reference['Depth (m)'], reference.theta, lower, upper)
        for date, group in selected[selected['Hole ID'] == well].groupby('Date'):
            group = group.sort_values('Depth (m)')
            reason = ''
            try:
                integral = integrate_profile(group['Depth (m)'], group.theta, lower, upper)
                col = f'well-{well}'
                depth = -gw.at[date,col] if col in gw and date in gw.index else np.nan
                if np.isfinite(depth) and depth >= 0.5 and depth <= upper:
                    raise ValueError('Water table within integration interval at survey date')
            except ValueError as exc:
                reason = str(exc)
            coverage.append(dict(well=well, date=date, lower_m=lower, upper_m=upper,
                                 measured_min_m=group['Depth (m)'].min(), measured_max_m=group['Depth (m)'].max(),
                                 included=not bool(reason), reason=reason))
            if not reason:
                records.append({'Hole ID':str(well), 'Date':date, 'Total Storage Change (mm)':integral-ref_integral})
    coverage = pd.DataFrame(coverage)
    coverage.to_csv(base/'neutron_probe_profile_qc.csv', index=False)
    domains = pd.DataFrame([dict(well=w, lower_m=a, upper_m=b) for w,(a,b) in DOMAINS.items()])
    counts = coverage.groupby('well').included.agg(n_surveys='size', n_retained='sum').reset_index()
    domains.merge(counts, on='well').to_csv(base/'neutron_probe_integration_domains.csv', index=False)
    storage = pd.DataFrame(records).sort_values(['Hole ID','Date'])
    storage['time_diff'] = storage.groupby('Hole ID').Date.diff()
    storage['is_new_segment'] = storage.time_diff > pd.Timedelta(days=60)
    storage['segment_id'] = storage.groupby('Hole ID').is_new_segment.cumsum()
    storage['line_group'] = storage['Hole ID'] + '_' + storage.segment_id.astype(str)
    storage.to_csv(base/'storage_ts.csv', index=False)
    print(f'Calibrated {len(cal)} observations; retained {len(storage)}/{len(coverage)} fixed-domain profiles.')
    print(domains.merge(counts,on='well').to_string(index=False))
    return storage


if __name__ == '__main__':
    process_neutron_profiles()
