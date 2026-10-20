# HBV sensitivity experiments

This self-contained module reproduces the rainfall-only HBV sensitivity experiments in Extended Data Fig. 5 of *Scale-dependent changes in global flood amplitude and concentration*.

From this directory:

```text
python -m pip install -r requirements.txt
python run_all.py
```

For the numerical simulation and FAI/FCI calculations without figure generation, only NumPy is needed:

```text
python run_all.py --skip-figures
```

The default output directory is `outputs/` beside these scripts. Use `--output-dir PATH` to select another destination. Existing output files with matching names are replaced. Do not run Python with `-O`, because the scientific checks include assertions.

No external datasets, downloads, parent-folder Python modules, GPU, or system fonts are required. The figure uses built-in Helvetica instead of platform-specific Arial; this changes typography only.

See [HBV_README.md](../../docs/HBV_README.md) for the equations, parameters, dependencies, expected files, measured run information and limitations.
