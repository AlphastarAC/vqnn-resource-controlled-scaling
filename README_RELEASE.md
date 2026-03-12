# Release-ready artifacts (GitHub / Zenodo)
This folder is prepared so you can publish the reproducibility artifacts with minimal extra work.

## Recommended workflow
1. Create a GitHub repository (private during review if needed).
2. Upload the entire `IEEE_Access_Final_Submission_Package/` directory or the contents of this folder.
3. After acceptance, create a tagged release (e.g., v1.0.0).
4. Connect the GitHub repo to Zenodo and generate a DOI for the release.
5. Replace placeholders in the manuscript:
   - https://github.com/<YOUR-GITHUB-HANDLE>/vqnn-resource-controlled-scaling
   - https://doi.org/10.5281/zenodo.<YOUR-RECORD-ID>

## Contents
- `REPRODUCIBILITY.md`: how to rerun the simulations/figures.
- `CITATION.cff`: citation metadata for GitHub.
- `zenodo.json`: metadata for Zenodo release.
- `LICENSE`: suggested license (MIT).
