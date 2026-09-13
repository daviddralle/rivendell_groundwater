# Rivendell weathered-bedrock groundwater analysis

Data and code for **Seasonal Activation and Drainage of a Weathered-Bedrock Groundwater System**.

The three numbered notebooks reproduce the manuscript analyses. Version 0.2.0 corrects neutron-probe calibration matching, uses fixed storage integration intervals, and restricts groundwater interpolation to complete short gaps. The manuscript observation window ends in WY2025.

## Reproduce the paper

Create the environment with `conda env create -f environment.yml` and activate `rivendell-groundwater`, or install `requirements.txt` in a Python environment. Run from this directory, in order:

1. `01_bedrock_vadose_storage_capacity.ipynb`: daily precipitation, neutron-probe calibration, fixed-domain storage, and storage-capacity fits.
2. `02_site_characterization.ipynb`: wet/dry comparison and annual precipitation–low-flow relationship.
3. `03_groundwater_streamflow_analysis.ipynb`: seasonal groundwater storage, weekly apparent coefficients, wet-season events, and seasonal groundwater peaks.

Run `python -m unittest -v test_analysis_utils.py` to check serial-specific calibration matching, fixed-domain integration, and interpolation of complete short gaps. The notebooks write data to `Data/` and figures to `Plots/`.

## Inputs and accounting

- `Data/raw/neutron_probe_rivendell_through_wy2025.csv` contains the Rivendell observations used by the paper, selected from the July 13, 2026 frozen export. The subset excludes other sites and observations after September 30, 2025. `source_row` retains the original export row number.
- `Data/raw/rivendell_borehole_stickup.csv` supplies casing geometry for the study wells. `Data/raw/probe_conversion_intervals.csv` supplies instrument serial numbers and calibration periods. `Data/raw/manuscript_source_manifest.json` documents selection and original-file hashes.
- `neutron_probe_processing.py` matches both serial number and calibration period, estimates inter-probe conversions from invariant depths, and applies casing-specific barrel calibrations. Repeated calibrated readings at a depth/date are combined by their median. No uncalibrated counts are silently accepted.
- Storage change is integrated relative to October 25, 2015 over fixed well-specific intervals, starting at 0.75 m beneath the reported soil profile. Incomplete profiles and internal gaps greater than 0.61 m are excluded; integration never extrapolates. `Data/neutron_probe_integration_domains.csv` gives the exact bounds and retained counts. These are monitored-interval estimates, not the entire vadose-zone storage capacity.
- Calibration coefficients, calibrated observations, rejected observations, and profile-coverage checks are exported as `Data/neutron_probe_*.csv`. The original broad source exports remain in the local workspace; the manuscript release contains only the selected source files above.
- `Data/angelo_precip_5m.csv` and `Data/prism_ppt.csv` generate the provenance-preserving `Data/site_precip_daily.csv`. Angelo days require at least 274 of 288 expected observations; PRISM supplies earlier years and incomplete days.
- `Data/groundwater_hourly_resampled.csv` and `Data/elder_creek_clean.csv` supply groundwater and daily runoff. Only complete interior groundwater gaps no longer than seven days are interpolated. `analysis_utils.py` implements the shared accounting rules.
- `Data/elder.tif` and `Data/wells.gpkg` supply stream-distance geometry and surveyed well elevations. Survey Z values embedded in the well layer are used for groundwater elevations; DEM elevations are retained only as an audit column.

## Main outputs

`storage_ts.csv`, `rivendell_theta_depth_timeseries.csv`, `vadose_storage_capacities.csv`, `seasonal_storage_balance.csv`, `groundwater_profile_geometry.csv`, `weekly_inferred_porosity.csv`, `recession_event_mass_balance.csv`, `low_flow_vs_water_year_precip.csv`, and `peak_groundwater_vs_precip.csv` provide the numerical results underlying the figures. The weekly file retains its legacy filename but contains **apparent storage coefficients**, not direct porosity measurements.

## Scope

The local checkout also contains earlier exploratory recession and matched-Q analyses (`04`, `05`, `06`, and their notes). They are preserved separately and are not required by, or included in, this manuscript release. Their existing outputs have not been regenerated for v0.2.0.

## Citation

Dralle, D. N., Rempe, D. M., Hahm, W. J., & Dietrich, W. E. (2026). *Data and code for: Seasonal Activation and Drainage of a Weathered-Bedrock Groundwater System* (v0.2.0). Zenodo. [All versions](https://doi.org/10.5281/zenodo.20361535).
