"""Reproduce the interannual recession analysis and manuscript figure.

Run after notebooks 01--03: python 07_interannual_recession.py
Uses the common six-well transect, one fixed reference, and observed response
criteria independent of streamflow. Fits ordinary least-squares regressions
and a pooled daily storage-flow relation.
"""
from pathlib import Path
import itertools
import json
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
DATA = OUT = ROOT/'Data'
PLOTS = ROOT/'Plots'
PLOTS.mkdir(exist_ok=True)
WELLS = ['well-1','well-3','well-5','well-7','well-10','well-15']
raw = pd.read_csv(DATA/'groundwater_hourly_resampled.csv',index_col=0,parse_dates=True)
gw = raw.resample('D').mean().where(raw.resample('D').count() >= 12)
rain = pd.read_csv(DATA/'site_precip_daily.csv',index_col=0,parse_dates=True).precip_mm
q = pd.read_csv(DATA/'elder_creek_clean.csv',index_col=0,parse_dates=True).flow_mm_day
if q.index.tz is not None:
    q.index = q.index.tz_localize(None)
d = gw.join(rain).join(q.rename('q'))
reference = gw.loc['2015-08-01':'2015-09-30',WELLS].median()
geom = pd.read_csv(DATA/'groundwater_profile_geometry.csv').set_index('well_id').loc[WELLS]
x = geom.dist_to_stream_m.to_numpy()
L = geom.loc['well-15','dist_to_stream_m']
weights = np.array([-2*L*L/3,1]) @ np.linalg.pinv(np.c_[x*x-2*L*x,np.ones(len(x))])
assert np.isclose(weights.sum(),1) and (weights > 0).all() and reference.notna().all()
d['S'] = (gw[WELLS]-reference) @ weights * 0.05 * 1000
pd.DataFrame({'well':WELLS,'distance_m':x,'weight':weights,
    'reference_signed_depth_m':reference.to_numpy()}).to_csv(
    OUT/'recession_reference_geometry.csv',index=False)
annual_P = pd.read_csv(DATA/'low_flow_vs_water_year_precip.csv').set_index('Water_Year')

def total(series,start,end):
    v = series.reindex(pd.date_range(start,end))
    return v.sum() if len(v) and v.notna().all() else np.nan

def state(t):
    v = d.reindex(pd.date_range(t-pd.Timedelta(days=1),t+pd.Timedelta(days=1)))
    return v.median().where(v.count() >= 2)

def events(year,start_month=1,end_month=6):
    first = pd.Timestamp(year,start_month,1)
    last = pd.Timestamp(year,end_month,1)+pd.offsets.MonthEnd(0)
    groups = []
    for t in rain.loc[first:last].index[rain.loc[first:last] > 1]:
        if not groups or (t-groups[-1][-1]).days > 3:
            groups.append([t])
        else:
            groups[-1].append(t)
    records = []
    for group in groups:
        start,end = group[0],group[-1]
        mm = total(rain,start,end)
        if not np.isfinite(mm) or mm < 10:
            continue
        before = d.reindex(pd.date_range(start-pd.Timedelta(days=3),start-pd.Timedelta(days=1)))
        baseline = before[WELLS].median().where(before[WELLS].count() >= 2)
        after = d.loc[start:end+pd.Timedelta(days=7)]
        rise = after[WELLS].max().where(after[WELLS].count() >= min(5,len(after)))-baseline
        s = d.loc[end+pd.Timedelta(days=1):end+pd.Timedelta(days=7),'S'].dropna()
        row = dict(WY=year,event_start=start,event_end=end,precip_mm=mm,
            available_wells=int(rise.notna().sum()),
            peak_date=s.idxmax() if len(s) >= 3 else pd.NaT)
        for threshold in [.03,.05,.10]:
            row[f'rises_{threshold:.2f}'] = int((rise >= threshold).sum())
        records.append(row)
    return records

all_events = pd.DataFrame([row for year in range(2015,2026) for row in events(year)])
all_events.to_csv(OUT/'recession_spring_response_events.csv',index=False)
rows,selection_rows = [],[]
for year in range(2015,2026):
    starts = {}
    yearly = all_events.loc[all_events.WY == year]
    for threshold in [.03,.05,.10]:
        selected = yearly.loc[(yearly.available_wells >= 4) & (yearly[f'rises_{threshold:.2f}'] >= 3)]
        start = selected.sort_values('event_end').iloc[-1].peak_date if len(selected) else pd.NaT
        if pd.notna(start) and ((yearly.event_start > start) & (yearly.available_wells < 4)).any():
            start = pd.NaT
        starts[threshold] = start
        selection_rows.append(dict(WY=year,rise_threshold_m=threshold,start_date=start))
    start = starts[.05]
    if pd.isna(start):
        continue
    st = state(start)
    sep = pd.Timestamp(year,9,1)
    end = state(sep)
    row = dict(WY=year,start_date=start,annual_P_mm=annual_P.loc[year,'Water_Year_Precip_mm'],
        S_start_mm=st.S,days_to_September=(sep-start).days,
        runoff_to_September_mm=total(q,start,sep),Q_September_mm_day=end.q,
        S_September_mm=end.S,P_start_to_September_mm=total(rain,start,sep))
    if not np.isfinite([row[k] for k in row if k not in ['WY','start_date']]).all():
        continue
    assert starts[.03] == start and starts[.10] == start
    rows.append(row)
pd.DataFrame(selection_rows).to_csv(OUT/'recession_start_selection.csv',index=False)
a = pd.DataFrame(rows).set_index('WY')
assert list(a.index) == [2015,2016,2017,2018,2019,2021,2022,2025]
summer = pd.DataFrame([row for year in a.index for row in events(year,7,8)])
assert not ((summer.available_wells >= 4) & (summer['rises_0.05'] >= 3)).any()
summer.to_csv(OUT/'recession_summer_response_audit.csv',index=False)

# Seasonal peak of complete-transect storage (November--May), for comparison
# with storage after the final detected response.
peak_rows = {}
for year in a.index:
    s = d.loc[f'{year-1}-11-01':f'{year}-05-31','S']
    s3 = s.rolling(3,center=True,min_periods=2).median()
    peak_rows[year] = dict(S_peak_mm=float(s3.max()),peak_date=s3.idxmax(),
        peak_complete_days=int(s.notna().sum()))
a = a.join(pd.DataFrame.from_dict(peak_rows,orient='index'))
a['storage_loss_mm'] = a.S_start_mm-a.S_September_mm
a['runoff_per_storage_loss'] = a.runoff_to_September_mm/a.storage_loss_mm

# Daily storage and flow from each selected start through September 1.
daily = pd.concat([d.loc[r.start_date:pd.Timestamp(y,9,1),['S','q']].dropna().assign(WY=y)
    for y,r in a.iterrows()])
daily['log_q'] = np.log(daily.q)
daily_fit = sm.OLS(daily.log_q,sm.add_constant(daily.S)).fit()
# Leave-one-year-out: predict each September 1 flow from that day's storage,
# using the daily relation fitted to the other seven recessions.
held = {}
for y in a.index:
    f = sm.OLS(daily.log_q[daily.WY != y],sm.add_constant(daily.S[daily.WY != y])).fit()
    held[y] = f.params['const']+f.params['S']*a.loc[y,'S_September_mm']
a['held_out_log_Q_September'] = pd.Series(held)
log_q_sep = np.log(a.Q_September_mm_day)
held_r2 = 1-((log_q_sep-a.held_out_log_Q_September)**2).sum()/((log_q_sep-log_q_sep.mean())**2).sum()
daily.rename_axis('date').to_csv(OUT/'recession_daily_storage_flow.csv')

specs = {
    'Precipitation to peak storage': ('S_peak_mm', ['annual_P_mm'], False),
    'Precipitation to initial storage': ('S_start_mm', ['annual_P_mm'], False),
    'Initial storage to runoff volume': ('runoff_to_September_mm', ['S_start_mm'], False),
    'Initial storage and later rainfall to runoff volume':
        ('runoff_to_September_mm', ['S_start_mm','P_start_to_September_mm'], False),
    'Initial storage and duration to low flow':
        ('Q_September_mm_day', ['S_start_mm','days_to_September'], True),
    'September storage to low flow': ('Q_September_mm_day', ['S_September_mm'], True),
    'Initial storage and duration to September storage':
        ('S_September_mm', ['S_start_mm','days_to_September'], False),
}

def permutation_p(y,X):
    """Exact permutation p value for R2: all orderings of the response."""
    X = sm.add_constant(X).to_numpy(float); y = np.asarray(y,float)
    H = X @ np.linalg.pinv(X)
    def r2(v):
        r = v-H@v
        return 1-r@r/((v-v.mean())@(v-v.mean()))
    obs = r2(y)
    perm = np.array([r2(y[list(p)]) for p in itertools.permutations(range(len(y)))])
    return float(np.mean(perm >= obs-1e-12))

fits, statistics, coefficients = {}, [], []
for name, (target, predictors, log_y) in specs.items():
    y = np.log(a[target]) if log_y else a[target]
    fit = sm.OLS(y, sm.add_constant(a[predictors])).fit()
    fits[name] = fit
    statistics.append(dict(model=name, response=target,
        response_transform='natural log' if log_y else 'none',
        n=int(fit.nobs), residual_df=int(fit.df_resid), R2=fit.rsquared,
        adjusted_R2=fit.rsquared_adj, F=fit.fvalue, model_p=fit.f_pvalue,
        permutation_p=permutation_p(y,a[predictors]),
        RMSE=np.mean(fit.resid**2)**.5, residual_SE=np.sqrt(fit.mse_resid)))
    intervals = fit.conf_int()
    for term in fit.params.index:
        coefficients.append(dict(model=name, term=term,
            estimate=fit.params[term], SE=fit.bse[term],
            t=fit.tvalues[term], p=fit.pvalues[term],
            CI95_lower=intervals.loc[term,0], CI95_upper=intervals.loc[term,1]))
# Daily values are serially correlated, so no p value is reported for this fit.
statistics.append(dict(model='Daily storage to daily flow, eight recessions pooled',
    response='q', response_transform='natural log', n=int(daily_fit.nobs),
    residual_df=int(daily_fit.df_resid), R2=daily_fit.rsquared,
    adjusted_R2=daily_fit.rsquared_adj, RMSE=np.mean(daily_fit.resid**2)**.5,
    residual_SE=np.sqrt(daily_fit.mse_resid)))
statistics.append(dict(model='Held-out September flow from daily relation (leave one year out)',
    response='Q_September_mm_day', response_transform='natural log', n=len(a), R2=held_r2,
    RMSE=float(np.mean((log_q_sep-a.held_out_log_Q_September)**2)**.5)))
for term in daily_fit.params.index:
    coefficients.append(dict(model='Daily storage to daily flow, eight recessions pooled',
        term=term, estimate=daily_fit.params[term]))
peak_fit = fits['Precipitation to peak storage']
storage_fit = fits['Precipitation to initial storage']
volume_fit = fits['Initial storage to runoff volume']
low_fit = fits['Initial storage and duration to low flow']
a['fitted_log_Q_September'] = low_fit.fittedvalues
# Direct exponentiation gives the fitted log-space central value, not a
# residual-smearing correction for an arithmetic conditional mean.
a['fitted_Q_September_mm_day'] = np.exp(low_fit.fittedvalues)
a['log_Q_residual'] = low_fit.resid
a.to_csv(OUT/'recession_reporting_values.csv')
pd.DataFrame(statistics).to_csv(OUT/'recession_reporting_statistics.csv', index=False)
pd.DataFrame(coefficients).to_csv(OUT/'recession_reporting_coefficients.csv', index=False)
beta = low_fit.params

plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':9,
    'axes.spines.top':False, 'axes.spines.right':False,
    'savefig.facecolor':'white', 'pdf.fonttype':42})

def ptext(p):
    return '$p < 0.001$' if p < .001 else f'$p = {p:.3f}$' if p < .1 else f'$p = {p:.2f}$'

# Manuscript figure: four panels at page width.
color = '#168b98'
year_colors = dict(zip(a.index,plt.get_cmap('tab10').colors))
fig,axs = plt.subplots(2,2,figsize=(7.4,6.9))
fig.subplots_adjust(left=.10,right=.98,bottom=.08,top=.95,wspace=.36,hspace=.42)
axs = axs.ravel()
titles = ['a) Precipitation and storage','b) Initial storage and runoff',
          'c) Daily storage and runoff','d) September 1 flow']
for ax,title in zip(axs,titles):
    ax.set_axisbelow(True); ax.grid(axis='y',color='.92',lw=.6)
    ax.tick_params(labelsize=8); ax.set_title(title,loc='left',fontsize=9.5,pad=8)

ax = axs[0]
ax.vlines(a.annual_P_mm,a.S_start_mm,a.S_peak_mm,color='.8',lw=.8,zorder=1)
ax.scatter(a.annual_P_mm,a.S_peak_mm,s=38,facecolors='white',edgecolors=color,linewidth=1.3,zorder=3,
    label=f'Seasonal peak: $R^2 = {peak_fit.rsquared:.2f}$; {ptext(peak_fit.f_pvalue)}')
ax.scatter(a.annual_P_mm,a.S_start_mm,s=38,color=color,edgecolors='#303030',linewidth=.5,zorder=4,
    label=f'After final response, $S_0$: $R^2 = {storage_fit.rsquared:.2f}$; {ptext(storage_fit.f_pvalue)}')
xx = np.linspace(a.annual_P_mm.min(),a.annual_P_mm.max(),50)
for fit,ls in [(peak_fit,'-'),(storage_fit,'--')]:
    ax.plot(xx,fit.params['const']+fit.params['annual_P_mm']*xx,ls,lw=1.1,color='.4',zorder=2)
off = {2015:(4,-11),2016:(4,-11),2017:(-27,-4),2018:(5,1),2019:(4,-11),2021:(4,-11),2022:(5,-10),2025:(4,-11)}
for year,row in a.iterrows():
    ax.annotate(str(year),(row.annual_P_mm,row.S_start_mm),xytext=off[year],textcoords='offset points',fontsize=7)
ax.set(xlim=(1000,3650),ylim=(100,390),xlabel='Water-year precipitation (mm)',ylabel='Storage anomaly (mm)')
ax.legend(loc='upper left',fontsize=7.2,frameon=False,handletextpad=.3,borderaxespad=.2)

ax = axs[1]
ax.scatter(a.S_start_mm,a.runoff_to_September_mm,s=38,color=color,edgecolors='#303030',linewidth=.5,zorder=3)
xx = np.linspace(a.S_start_mm.min(),a.S_start_mm.max(),50)
ax.plot(xx,volume_fit.params['const']+volume_fit.params['S_start_mm']*xx,'--',lw=1.1,color='.4',zorder=1)
off = {2015:(-30,3),2016:(4,-12),2017:(-31,6),2018:(5,-12),2019:(4,-13),2021:(4,7),2022:(4,-13),2025:(-30,3)}
for year,row in a.iterrows():
    ax.annotate(str(year),(row.S_start_mm,row.runoff_to_September_mm),xytext=off[year],textcoords='offset points',fontsize=7)
ax.text(.03,.97,f'$R^2 = {volume_fit.rsquared:.2f}$; {ptext(volume_fit.f_pvalue)}\n'
    f'slope = {volume_fit.params["S_start_mm"]:.2f} mm mm$^{{-1}}$',transform=ax.transAxes,va='top',fontsize=8,linespacing=1.6)
ax.set(xlim=(115,250),ylim=(25,270),xlabel='Initial storage anomaly, $S_0$ (mm)',ylabel='Runoff to September 1 (mm)')

ax = axs[2]
for year,g in daily.groupby('WY'):
    ax.plot(g.S,g.q,lw=1.2,color=year_colors[year],label=str(year),zorder=2)
ax.scatter(a.S_September_mm,a.Q_September_mm_day,s=22,c=[year_colors[y] for y in a.index],
    edgecolors='#303030',linewidth=.6,zorder=4)
xx = np.linspace(daily.S.min(),daily.S.max(),50)
ax.plot(xx,np.exp(daily_fit.params['const']+daily_fit.params['S']*xx),'--',lw=1.1,color='.25',zorder=3)
ax.set_yscale('log')
ax.text(.97,.04,f'{len(daily)} days; $R^2 = {daily_fit.rsquared:.2f}$ (log flow)',transform=ax.transAxes,ha='right',fontsize=8)
ax.set(xlabel='Storage anomaly (mm)',ylabel='Runoff (mm/day)')
ax.legend(loc='upper left',fontsize=7,frameon=False,ncol=2,handlelength=1.2,columnspacing=.8,borderaxespad=.2)

ax = axs[3]
ax.scatter(a.Q_September_mm_day,a.fitted_Q_September_mm_day,s=38,color=color,edgecolors='#303030',linewidth=.5,zorder=3)
off = {2015:(-28,-13),2016:(4,-13),2017:(4,5),2018:(5,-13),2019:(-28,5),2021:(4,5),2022:(-30,0),2025:(4,5)}
for year,row in a.iterrows():
    ax.annotate(str(year),(row.Q_September_mm_day,row.fitted_Q_September_mm_day),xytext=off[year],textcoords='offset points',fontsize=7)
ax.plot([.04,.205],[.04,.205],'--',color='.4',lw=1.1,zorder=1)
ax.text(.03,.97,f'Initial storage and duration\n$R^2 = {low_fit.rsquared:.2f}$ (log flow); {ptext(low_fit.f_pvalue)}',
    transform=ax.transAxes,va='top',fontsize=8,linespacing=1.6)
ax.set(xlim=(.04,.205),ylim=(.04,.205),xlabel='Observed September 1 runoff (mm/day)',ylabel='Fitted September 1 runoff (mm/day)')
ax.set_xticks([.05,.10,.15,.20]); ax.set_yticks([.05,.10,.15,.20])
for ax in axs:
    ax.xaxis.label.set_size(9); ax.yaxis.label.set_size(9)
for ext in ['pdf','png']:
    fig.savefig(PLOTS/f'recession_reporting_manuscript.{ext}',dpi=250)
plt.close(fig)

metadata = {
    'fit':'Ordinary least squares, intercept included; all eight years used for the annual fits.',
    'permutation_p':'Exact: share of all 40320 orderings of the response with R2 at least as large as observed.',
    'low_flow_response':'Natural logarithm of September 1 flow (mm/day).',
    'low_flow_display':'Observed versus exp(fitted log flow); R2 is computed in log space.',
    'low_flow_model':'ln(Q/(1 mm/day)) = a + b*S0_mm + c*T_days',
    'low_flow_coefficients':{k:float(v) for k,v in beta.items()},
    'initial_storage':'Three-day median centered on the selected post-event maximum; fixed six-well transect and 0.05 effective drainable porosity; one common 2015 datum.',
    'peak_storage':'Maximum centered three-day median of complete six-well storage, November 1 through May 31; peak_complete_days counts days with all six wells.',
    'storage_loss':'Initial storage minus the September 1 three-day median storage.',
    'daily_relation':'Natural-log daily flow against daily storage from each selected start through September 1, pooled over eight years. Daily values are serially correlated; the fit is descriptive. The held-out test predicts each September 1 flow from that day\'s storage using the other seven years.',
    'start':'Final detected groundwater response; see manuscript Methods and SI Text S3. This does not declare subsequent recharge zero.',
    'duration':'Calendar days from selected start through September 1; final-response date and duration are not separate predictors.',
    'runoff_volume':'Daily runoff summed inclusively from selected start through September 1.',
    'September_flow':'Median daily flow over August 31-September 2.',
    'precipitation':'Water-year total; a retrospective descriptor.',
    'interpretation':'Hillslope storage and catchment runoff; fitted associations do not imply strict control-volume closure.',
    'years':[int(v) for v in a.index],
}
(OUT/'recession_reporting_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
pd.set_option('display.width',250)
print(pd.DataFrame(statistics).round(4).to_string(index=False))
print(pd.DataFrame(coefficients).round(5).to_string(index=False))
print(a[['S_peak_mm','peak_date','peak_complete_days','S_start_mm','S_September_mm','storage_loss_mm','runoff_per_storage_loss','held_out_log_Q_September']].round(3).to_string())

# Supporting table is generated from the exact fitted values.
table = r"""\begin{table}[htbp]
\centering
\caption{Annual values used in the interannual recession analysis. $P$ is water-year precipitation; $S_{peak}$ is the seasonal peak storage anomaly; $S_0$ is the storage anomaly after the final detected response; $T$ is days from that start to September 1; $S_{Sep}$ is the September 1 storage anomaly; $V$ is cumulative runoff over the interval; $Q_{Sep}$ is September 1 flow; and $P_{after}$ is precipitation from the start through September 1. Storage uses a common fixed reference.}
\label{tab:annual_recession}
\small
\setlength{\tabcolsep}{4pt}
\begin{tabular}{rlrrrrrrrr}
\hline
WY & Start & $P$ & $S_{peak}$ & $S_0$ & $T$ & $S_{Sep}$ & $V$ & $Q_{Sep}$ & $P_{after}$ \\
 & & mm & mm & mm & days & mm & mm & mm d$^{-1}$ & mm \\
\hline
"""
for year,row in a.iterrows():
    table += f'{year} & {row.start_date.strftime("%b %d")} & {row.annual_P_mm:.0f} & {row.S_peak_mm:.1f} & {row.S_start_mm:.1f} & {row.days_to_September:.0f} & {round(row.S_September_mm,1)+0:.1f} & {row.runoff_to_September_mm:.1f} & {row.Q_September_mm_day:.3f} & {row.P_start_to_September_mm:.1f}'+r' \\'+'\n'
table += r'\hline'+'\n'+r'\end{tabular}'+'\n'+r'\end{table}'+'\n'
(OUT/'recession_annual_table.tex').write_text(table)
