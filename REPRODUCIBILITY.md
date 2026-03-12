# Reproducibility (offline + online)
This submission includes an offline supplementary package that can reproduce all results without internet access.

## Offline (during review)
1. Unzip the submission package.
2. Install dependencies listed in `requirements.txt` (Python >= 3.10 recommended).
3. Run `python reproduce_all.py` to regenerate:
   - Tables (grid + coefficients + budget study)
   - Figures (Fig. 1–6)
4. Confirm the regenerated files match `figs/` and `supplement/` outputs.

## Online (after acceptance)
1. Push this repository to GitHub.
2. Tag a release (e.g., v1.0.0).
3. Create a Zenodo archive for that release and obtain a DOI.
4. Update the manuscript Data/Code Availability section with the final URLs.
