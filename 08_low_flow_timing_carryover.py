"""Timing and carryover checks on the 35-year precipitation and runoff record.

Run after notebook 02: python 08_low_flow_timing_carryover.py
Tests whether September 1 flow and annual minimum runoff depend on the previous
year's precipitation, on April--June precipitation, and on days since the last
large storm. Uses no groundwater data. Ordinary least squares on natural logs.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parent
DATA = ROOT/'Data'
rain = pd.read_csv(DATA/'site_precip_daily.csv',index_col=0,parse_dates=True).precip_mm
q = pd.read_csv(DATA/'elder_creek_clean.csv',index_col=0,parse_dates=True).flow_mm_day
if q.index.tz is not None:
    q.index = q.index.tz_localize(None)
df = pd.read_csv(DATA/'low_flow_vs_water_year_precip.csv').set_index('Water_Year')

def last_storm_end(year,threshold):
    """End of the last January--August event totaling at least `threshold` mm.
    Wet days exceed 1 mm; wet days no more than three days apart are grouped."""
    r = rain.loc[pd.Timestamp(year,1,1):pd.Timestamp(year,8,31)]
    groups = []
    for t in r.index[r > 1]:
        if not groups or (t-groups[-1][-1]).days > 3:
            groups.append([t])
        else:
            groups[-1].append(t)
    ends = [g[-1] for g in groups if rain.loc[g[0]:g[-1]].sum() >= threshold]
    return ends[-1] if ends else pd.NaT

df['ln_P'] = np.log(df.Water_Year_Precip_mm)
df['ln_P_previous'] = df.ln_P.shift(1)
df['ln_Q_min'] = np.log(df.Minimum_Runoff_mm_day)
df['Q_September_mm_day'] = [q.reindex(pd.date_range(f'{y}-08-31',f'{y}-09-02')).median() for y in df.index]
df['ln_Q_September'] = np.log(df.Q_September_mm_day)
df['P_April_June_mm'] = [rain.loc[f'{y}-04-01':f'{y}-06-30'].sum() for y in df.index]
df['ln_P_April_June'] = np.log(df.P_April_June_mm)
THRESHOLDS = [10,25,50,75,100]
for thr in THRESHOLDS:
    df[f'days_since_storm_{thr}mm'] = [(pd.Timestamp(y,9,1)-last_storm_end(y,thr)).days for y in df.index]

rows = []
def fit(response,predictors):
    m = sm.OLS(df[response],sm.add_constant(df[predictors]),missing='drop').fit()
    for term in predictors:
        lo,hi = m.conf_int().loc[term]
        rows.append(dict(response=response,predictors=' + '.join(predictors),n=int(m.nobs),
            R2=m.rsquared,adjusted_R2=m.rsquared_adj,term=term,estimate=m.params[term],
            CI95_lower=lo,CI95_upper=hi,p=m.pvalues[term]))

for response in ['ln_Q_September','ln_Q_min']:
    fit(response,['ln_P'])
    fit(response,['ln_P','ln_P_previous'])
    fit(response,['ln_P_April_June'])
    fit(response,['ln_P','ln_P_April_June'])
    for thr in THRESHOLDS:
        fit(response,[f'days_since_storm_{thr}mm'])
        fit(response,['ln_P',f'days_since_storm_{thr}mm'])
out = pd.DataFrame(rows)
for thr in THRESHOLDS:
    t = df[f'days_since_storm_{thr}mm']
    out = pd.concat([out,pd.DataFrame([dict(response=f'days_since_storm_{thr}mm',predictors='ln_P (Pearson r)',
        n=int(t.notna().sum()),term='r',estimate=df.ln_P.corr(t))])])
out.to_csv(DATA/'low_flow_timing_carryover.csv',index=False)
df.to_csv(DATA/'low_flow_timing_carryover_values.csv')
pd.set_option('display.width',250); pd.set_option('display.max_rows',200)
print(out.round(4).to_string(index=False))
print('April-June precipitation range (mm):',df.P_April_June_mm.min().round(1),df.P_April_June_mm.max().round(1))
