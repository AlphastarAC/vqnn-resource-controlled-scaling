# Supplementary package (IEEE Access submission)

This folder contains the rendered figures and CSV summaries used to build the manuscript.

## Contents
- `main.tex`, `main.pdf`: manuscript source and compiled PDF.
- `figs/`: figure image files.
- `supplement/`: CSV tables referenced by the manuscript.
- `reproduce_monte_carlo.py`: fixed-seed script that reproduces the robustness post-processing artifacts:
  - `figs/fig8_mc_gain.png`
  - `supplement/mc_gain_summary.csv`
  - `supplement/mc_coefficients_summary.csv`
  - `figs/fig9_p_calibration_gain.png`
  - `supplement/p_calibration_gain_summary.csv`
  - `figs/fig10_budget_proxy.png`
  - `supplement/budget_proxy_optimization.csv`
  - `supplement/table_budget_proxy_p002.csv`

## Reproducing the Monte Carlo artifacts
From the package root:

```bash
python reproduce_monte_carlo.py
```

Dependencies: Python 3, NumPy, Pandas, Matplotlib.

The script uses a fixed random seed to make the results deterministic.
