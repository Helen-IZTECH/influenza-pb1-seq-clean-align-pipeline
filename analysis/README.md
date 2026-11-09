# Influenza A PB1 / PB1–PB2 Interface — Conservation & Druggability Pipeline

A fully computational analysis of the influenza A **PB1–PB2 polymerase interface**
as a resistance-robust antiviral drug target. The pipeline goes from raw PB1
coding sequences to population-genetics statistics, per-residue conservation, and
a focused structural analysis of the PB1 C-terminus / PB2 N-terminus binding
pocket.

## Data

- **H1N1pdm09 PB1 coding sequences, 2009–2025, 5 continents** (GISAID / NCBI
  Virus), aligned to the **A/California/04/2009** reference (ACP41103.1) —
  gapless, in-frame, exactly 2274 nt.
- **PB1 sequences for additional subtypes** (H7N9, H9N2, H5N1, H7N7) for the
  cross-subtype comparison.
- **Structure: PDB 3A1G** (PB1 C-terminus / PB2 N-terminus interface).

> Sequence data and structures are **not** included in this repo. Point each
> script's `CONFIG` block at your local copies.

## Repository layout

The pipeline runs in three ordered stages.

### `R/01_data_pipeline/` — from raw downloads to a clean dataset
| Script | What it does |
|---|---|
| `01_standardize_names.R` | Walks the raw dataset and copies every FASTA into a standardized `PB1_{YEAR}/H1N1_PB1_{CONTINENT}_{YEAR}` tree. Originals are never modified. |
| `02_filter_perfect_cds.R` | Keeps only "perfect" CDS (2274 nt, ATG start, valid stop, no premature stop, no N); writes rejects and a filtering report. |
| `03_qc_report.R` | Per-continent QC: gaps, N content, frameshift, premature stop → styled Excel report. |

### `R/02_popgen_conservation/` — conservation & population genetics
| Script | What it does |
|---|---|
| `01_site_by_site_conservation.R` | Per-position nucleotide (2274) and amino-acid (757) conservation: consensus, global %, Shannon entropy, per-year, year-to-year shift, per-continent. |
| `02_dnasp_engine.R` | DnaSP-equivalent engine: segregating sites, π, Watterson's θ, haplotype diversity, syn/nonsyn polymorphism, Tajima's D, Fu & Li's D*/F*. |
| `03_dnasp_boxplots.R` | Boxplots of the DnaSP metrics by year and by continent. |
| `04_master_workbook.R` | Assembles a 9-sheet master workbook, incl. the Sugiyama et al. 2009 PB1–PB2 contact-pocket sheet. |
| `05_domain_figures.R` | Four publication figures (domain dumbbell/lollipop, domain×year heatmap, PB1–PB2 interface + continents, within-domain boxplot) + multi-page PDF. |

### `R/03_interface_analysis/` — the PB1–PB2 interface story
| Script | What it does |
|---|---|
| `01_contact_conservation.R` | Defines the 5 Å contact set, per-residue entropy, invariant core vs evolvable residues. Figure 1. |
| `02_temporal_geographic.R` | Substitution over time and across continents (positions 757, 729). Figure 2. |
| `03_linkage_phylogeny.R` | Co-occurrence of 757R/729E, a year-matched test, NJ tree coloured by residue state. Figure 3. |
| `04_pp7_subtypes.R` | PP7 pocket (+8 Å neighbours) and whether it differs between PP7-sensitive and PP7-resistant subtypes. Figure 4. |

Within stage 3, scripts 02–04 reuse objects from `01_contact_conservation.R`, so
run them in order.

Two further stages live in sub-folders and are run independently (each has its own README):

| Folder | What it does |
|---|---|
| `05_selection_dnds/` | Site-wise dN/dS (HyPhy FEL / FUBAR / MEME) on a codon alignment; per-site CSVs. Setup via `environment.yml`, run with `run_selection.sh`. |
| `06_ddg_binding/` | EvoEF2 binding ΔΔG scan of the interface, for **PB1** (`bash run_ddg.sh PB1`) or **PB2** (`bash run_ddg.sh PB2`); joins the stage-05 selection tables. Works on Linux, WSL and macOS. |

Order for the selection + energetics layer: stage 05 (selection tables) can be run
before or after 06 (ΔΔG); they only meet in `06_ddg_binding/merge_ddg_selection.R`.

## Conventions

Each script starts with a `CONFIG` block — edit the input/output paths there.
The expected working layout is:

```
data/       raw and standardized sequence trees
results/    Excel tables (conservation, DnaSP, master workbook)
figures/    PNG / PDF figures
```

## Requirements

R (≥ 4.1) with: `Biostrings`, `bio3d`, `ape`, `pegas`, `phangorn`, `ggplot2`,
`dplyr`, `tidyr`, `patchwork`, `scales`, `gridExtra`, `readxl`, `openxlsx`,
`stringi`. Scripts install any missing CRAN packages on first run; `Biostrings`
is from Bioconductor (`BiocManager::install("Biostrings")`).

## Reference models

- Sugiyama et al. 2009 (PB1–PB2 interface contact residues)
- Pflug et al. 2014 (PB1 structural domain model)
- PDB 3A1G (interface structure)
