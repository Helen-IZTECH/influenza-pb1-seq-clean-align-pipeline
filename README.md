
# Influenza PB1 sequence preprocessing pipeline (R + MAFFT)

Small R scripts developed to support an undergraduate project on Influenza A PB1 sequence analysis.

## What this does
Given a raw FASTA file, the main script:
- reads sequences and performs basic QC (length, invalid characters, gap checks, N content)
- extracts a fixed-length CDS starting from the first ATG (default: 2274 nt)
- runs MAFFT for alignment and trims all-gap columns at the alignment ends
- writes:
  - `<prefix>_final.fasta` (aligned + trimmed)
  - `<prefix>_QC.csv` (QC summary table)
- performs a final “DnaSP-ready” cleaning step and writes:
  - `<prefix>_final_clean.fasta`
  - `<prefix>_removed_ids.txt`

## Requirements
- R
- Biostrings (Bioconductor)
- MAFFT installed locally (script calls MAFFT via `system2()`)

## How to run (example)
1. Open `scripts/pb1_pipeline.R`
2. Set:
   - working directory (where your FASTA is)
   - `input_fasta`
   - `output_prefix`
   - `mafft_path`
3. Run the script in R/RStudio.

## Country-balanced sampling utility

- scripts/pb1_country_balanced_sampling.R  
  Country-balanced round-robin subsampling from FASTA files to reduce geographic over-representation.


This script implements a simple round-robin sampling strategy to select a fixed number of sequences
from large PB1 FASTA files while limiting over-representation from any single country.

The goal is not optimization, but bias control for downstream evolutionary and population-genetic analyses.

Key features:
- country-aware round-robin selection
- adaptive per-country cap to reach target sample size
- reproducible sampling via fixed random seed
- single FASTA output plus in-R summary tables for sanity checks

This utility was developed to support exploratory analyses on large influenza PB1 datasets.

## Notes
- This repository focuses on workflow logic and reproducible outputs for research use,
  not on production-grade software engineering.
- Input data are not included.
