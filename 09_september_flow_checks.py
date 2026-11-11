"""Checks on the September 1 endpoint and on peak-to-start storage changes.

Run after 07 and 08: python 09_september_flow_checks.py
(1) How closely September 1 flow tracks the annual minimum over the 35-year
    record, and when the minimum occurs.
(2) Runoff between September 1 and the annual minimum in the eight recession
    years, as a depth and as a fraction of start-to-September runoff.
(3) Storage change and elapsed days between the seasonal peak and the final
    detected response, and storage loss minus runoff, for the same years.
Writes Data/september_flow_checks.csv and Data/september_flow_checks_annual.csv.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = ROOT/'Data'
q = pd.read_csv(DATA/'elder_creek_clean.csv',index_col=0,parse_dates=True).flow_mm_day
if q.index.tz is not None:
    q.index = q.index.tz_localize(None)
lf = pd.read_csv(DATA/'low_flow_timing_carryover_values.csv',index_col=0,parse_dates=['Minimum_Date'])
rec = pd.read_csv(DATA/'recession_reporting_values.csv',index_col=0,parse_dates=['start_date','peak_date'])

# (1) September 1 flow against the annual minimum, WY1991-2025
ratio = lf.Q_September_mm_day/lf.Minimum_Runoff_mm_day
doy = lf.Minimum_Date.dt.dayofyear
summary = dict(
    n_years=len(lf),
    r_log_QSep_vs_log_Qmin=np.corrcoef(lf.ln_Q_September,lf.ln_Q_min)[0,1],
    QSep_over_Qmin_median=ratio.median(),QSep_over_Qmin_min=ratio.min(),QSep_over_Qmin_max=ratio.max(),
    minimum_date_median=(pd.Timestamp(2001,1,1)+pd.Timedelta(days=int(doy.median())-1)).strftime('%b %d'),
    minimum_date_q25=(pd.Timestamp(2001,1,1)+pd.Timedelta(days=int(doy.quantile(.25))-1)).strftime('%b %d'),
    minimum_date_q75=(pd.Timestamp(2001,1,1)+pd.Timedelta(days=int(doy.quantile(.75))-1)).strftime('%b %d'),
    years_minimum_after_Sep1=int((doy > pd.Timestamp(2001,9,1).dayofyear).sum()))

# (2), (3) the eight recession years
rows = []
for wy,r in rec.iterrows():
    tmin = lf.Minimum_Date[wy]
    after = q.loc[pd.Timestamp(wy,9,2):tmin].sum() if tmin > pd.Timestamp(wy,9,1) else 0.0
    rows.append(dict(WY=wy,minimum_date=tmin.date(),runoff_Sep1_to_minimum_mm=after,
        fraction_of_start_to_Sep_runoff=after/r.runoff_to_September_mm,
        peak_to_start_storage_change_mm=r.S_peak_mm-r.S_start_mm,
        peak_to_start_days=(r.start_date-r.peak_date).days,
        storage_loss_minus_runoff_mm=r.storage_loss_mm-r.runoff_to_September_mm,
        excess_fraction_of_storage_loss=(r.storage_loss_mm-r.runoff_to_September_mm)/r.storage_loss_mm))
annual = pd.DataFrame(rows).set_index('WY')
summary.update(
    runoff_Sep1_to_minimum_mm_min=annual.runoff_Sep1_to_minimum_mm.min(),
    runoff_Sep1_to_minimum_mm_max=annual.runoff_Sep1_to_minimum_mm.max(),
    fraction_of_V_min=annual.fraction_of_start_to_Sep_runoff.min(),
    fraction_of_V_max=annual.fraction_of_start_to_Sep_runoff.max(),
    S_peak_slope_mm_per_100mm_P=100*np.polyfit(rec.annual_P_mm,rec.S_peak_mm,1)[0],
    P_range_fold=lf.Water_Year_Precip_mm.max()/lf.Water_Year_Precip_mm.min(),
    Qmin_range_fold=lf.Minimum_Runoff_mm_day.max()/lf.Minimum_Runoff_mm_day.min())
pd.Series(summary,name='value').rename_axis('quantity').to_csv(DATA/'september_flow_checks.csv')
annual.to_csv(DATA/'september_flow_checks_annual.csv')
pd.set_option('display.width',250)
print(pd.Series(summary).to_string()); print(annual.round(3).to_string())
