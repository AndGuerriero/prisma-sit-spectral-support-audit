# PRISMA-SIT spectral-support audit

This repository contains the MATLAB and Python code used for the manuscript:

**Field-Level Spectral Support Auditing of Legacy Land-Use Inventories Using PRISMA Hyperspectral Data**

The workflow audits legacy SIT land-use labels using PRISMA hyperspectral imagery. The method estimates class-specific spectral prototypes from the assigned original SIT LIVELLO_4 labels, computes pixel-level spectral-support scores, aggregates them at field-scene level, and assigns diagnostic spectral-support strata.

## Repository contents

- `matlab/`: MATLAB scripts used to extract PRISMA-SIT pixel samples and export scene-balanced datasets.
- `notebooks/`: Jupyter notebooks used for the main audit and sensitivity analyses.
- `scripts/`: Python scripts for reusable processing and plotting functions.
- `tables/`: LaTeX tables generated for the manuscript.
- `figures/`: figures generated from the analysis.
- `metadata/`: lightweight metadata files needed for documentation and plotting.

## Data availability

The PRISMA hyperspectral data are available from the Italian Space Agency (ASI) through the official PRISMA portal, subject to ASI access policies.

The original SIT inventory data are not redistributed in this repository. Users should obtain the source inventory from the official provider.

Large intermediate files, including per-scene spectral matrices, are not included because of size and data-access constraints.

## Main workflow

1. Export valid PRISMA-SIT pixels using the MATLAB scripts in `matlab/`.
2. Run the main spectral-support audit notebook.
3. Run sensitivity analyses:
   - SAD/EUC distance sensitivity
   - operating-threshold sensitivity
   - leave-one-scene-out prototype sensitivity
   - cloudiest-scene sensitivity
   - per-scene cap seed sensitivity
4. Generate manuscript tables and figures.

## Python environment

Install Python dependencies with:

```bash
pip install -r requirements.txt