Reproduce: run run_all.py from the HBV code directory (Python, NumPy, ReportLab, pypdfium2 and Pillow). See package/docs/HBV_README.md.
The CSV files contain all daily values, states/fluxes, scenario settings/metrics and routing weights.
validation.json documents full-record mass conservation and finite-record residual storage.
index_fit_validation.json documents the separate 1-20-day FAI/FCI fitting; fitted_indices.csv contains all indices and R-squared values.
maximum_window_volumes.csv records the maximizing windows over the full source series, independently of the displayed 0-20-day interval.
methods.txt specifies equations, flux ordering, parameters, assumptions and sources.
No observations, external datasets, random samples or calibration were used.
The resulting hydrographs are simulated experiments, not observed flood events.
