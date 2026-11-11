# Frozen manuscript inputs

The manuscript uses `neutron_probe_rivendell_through_wy2025.csv`, `rivendell_borehole_stickup.csv`, and `probe_conversion_intervals.csv`. These were selected from the July 13, 2026 source exports; the observation subset includes only study wells at Rivendell through September 30, 2025. `source_row` preserves row numbers in the original observation export, including its header. Selection rules are recorded in `manuscript_source_manifest.json`.

The full original exports remain unchanged in the local workspace but are not needed in the manuscript release. Regeneration reads the selected files directly and requires no Google Sheets connection.
