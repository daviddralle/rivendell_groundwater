# Rivendell weathered-bedrock groundwater analysis

Data and code for **From Storms to Low Flows: Storage and Drainage in a Bedrock Hillslope Aquifer**.

The three numbered notebooks and scripts `07`--`09` reproduce the manuscript analyses. Version 0.2.0 corrected neutron-probe calibration matching, used fixed storage integration intervals, and restricted groundwater interpolation to complete short gaps. Version 0.3.0 adds the interannual recession analysis. The manuscript observation window ends in WY2025.

## Reproduce the paper

Create the environment with `conda env create -f environment.yml` and activate `rivendell-groundwater`, or install `requirements.txt` in a Python environment. Run from this directory, in order:

1. `01_bedrock_vadose_storage_capacity.ipynb`: daily precipitation, neutron-probe calibration, fixed-domain storage, and storage-capacity fits.
2. `02_site_characterization.ipynb`: wet/dry comparison and the annual precipitation–low-flow relationship, including the precipitation elasticity of low flow.
3. `03_groundwater_streamflow_analysis.ipynb`: seasonal groundwater storage, weekly apparent coefficients, wet-season events, and seasonal groundwater peaks.
4. `python 07_interannual_recession.py`: response-based recession starts, seasonal peak and initial groundwater storage, subsequent runoff, the pooled daily storage–runoff relation, and the fitted September 1 low-flow model. It also writes the manuscript figure and SI annual-values table.
5. `python 08_low_flow_timing_carryover.py`: 35-year checks of September 1 flow and minimum runoff against previous-year precipitation, April–June precipitation, and days since the last large storm. Writes `Data/low_flow_timing_carryover.csv` and `Data/low_flow_timing_carryover_values.csv`.
6. `python 09_september_flow_checks.py`: how closely September 1 flow tracks the annual minimum over the 35-year record, runoff between September 1 and the minimum in the eight recession years, and the storage change and elapsed days between the seasonal peak and the final detected response. Writes `Data/september_flow_checks.csv` and `Data/september_flow_checks_annual.csv`.

Run `python -m unittest -v test_analysis_utils.py` to check serial-specific calibration matching, fixed-domain integration, and interpolation of complete short gaps. The notebooks write data to `Data/` and figures to `Plots/`.

## Inputs and accounting

- `Data/raw/neutron_probe_rivendell_through_wy2025.csv` contains the Rivendell observations used by the paper, selected from the July 13, 2026 frozen export. The subset excludes other sites and observations after September 30, 2025. `source_row` retains the original export row number.
- `Data/raw/rivendell_borehole_stickup.csv` supplies casing geometry for the study wells. `Data/raw/probe_conversion_intervals.csv` supplies instrument serial numbers and calibration periods. `Data/raw/manuscript_source_manifest.json` documents the selection.
- `neutron_probe_processing.py` matches both serial number and calibration period, estimates inter-probe conversions from invariant depths, and applies casing-specific barrel calibrations. Repeated calibrated readings at a depth/date are combined by their median. No uncalibrated counts are silently accepted.
- Storage change is integrated relative to October 25, 2015 over fixed well-specific intervals, starting at 0.75 m beneath the reported soil profile. Incomplete profiles and internal gaps greater than 0.61 m are excluded; integration never extrapolates. `Data/neutron_probe_integration_domains.csv` gives the exact bounds and retained counts. These are monitored-interval estimates, not the entire vadose-zone storage capacity.
- Calibration coefficients, calibrated observations, rejected observations, and profile-coverage checks are exported as `Data/neutron_probe_*.csv`. The original broad source exports remain in the local workspace; the manuscript release contains only the selected source files above.
- `Data/angelo_precip_5m.csv` and `Data/prism_ppt.csv` generate the provenance-preserving `Data/site_precip_daily.csv`. Angelo days require at least 274 of 288 expected observations; PRISM supplies earlier years and incomplete days.
- `Data/groundwater_hourly_resampled.csv` and `Data/elder_creek_clean.csv` supply groundwater and daily runoff. Only complete interior groundwater gaps no longer than seven days are interpolated. `analysis_utils.py` implements the shared accounting rules.
- `Data/elder.tif` and `Data/wells.gpkg` supply stream-distance geometry and surveyed well elevations. Survey Z values embedded in the well layer are used for groundwater elevations; DEM elevations are retained only as an audit column.
- The interannual analysis uses a fixed six-well transect, daily means requiring at least 12 hourly observations, and no further gap filling. Storage uses one common August--September 2015 groundwater reference. Its final-response criterion is independent of streamflow; retained years and event checks are written with the results. This groundwater reference is distinct from the October 25, 2015 neutron-probe reference.

## Main outputs

`storage_ts.csv`, `rivendell_theta_depth_timeseries.csv`, `vadose_storage_capacities.csv`, `seasonal_storage_balance.csv`, `groundwater_profile_geometry.csv`, `weekly_inferred_porosity.csv`, `recession_event_mass_balance.csv`, `low_flow_vs_water_year_precip.csv`, and `peak_groundwater_vs_precip.csv` provide the numerical results underlying the figures. The weekly file retains its legacy filename but contains **apparent storage coefficients**, not direct porosity measurements.

The interannual analysis writes `recession_reporting_values.csv`, `recession_reporting_statistics.csv`, `recession_reporting_coefficients.csv`, `recession_daily_storage_flow.csv`, and `recession_reporting_metadata.json`. `Plots/recession_reporting_manuscript.pdf` is the manuscript figure and `Data/recession_annual_table.tex` the SI table. Annual regressions are ordinary least-squares fits using all eight complete years, with exact permutation p values; the low-flow response is natural-log flow. The pooled daily fit is descriptive because daily values are serially correlated. `low_flow_precip_elasticity.csv` holds the 35-year correlation and elasticity.

## Scope

The local checkout also contains earlier exploratory recession and matched-Q analyses (`04`, `05`, `06`, and their notes). They are preserved separately and are not required by, or included in, this manuscript release. Their existing outputs have not been regenerated.

## Citation

Dralle, D. N., Rempe, D. M., Hahm, W. J., & Dietrich, W. E. (2026). *Data and code for: From Storms to Low Flows: Storage and Drainage in a Bedrock Hillslope Aquifer* (v0.3.0). Zenodo. [All versions](https://doi.org/10.5281/zenodo.20361535).
