# 05 — Site-wise dN/dS selection analysis (HyPhy)

Tests whether the PB1–PB2 interface is under purifying (negative) or positive
selection, site by site, using three independent methods on a 500-sequence
representative H1N1pdm09 set (2009–2025).

- **FEL** — pervasive selection, frequentist (per-site likelihood-ratio test)
- **FUBAR** — pervasive selection, Bayesian (per-site posterior probability)
- **MEME** — episodic positive selection (adaptive change on a subset of branches)

## Files

| File | Purpose |
|---|---|
| `environment.yml` | Conda recipe (hyphy, iqtree, python) — reproduces the exact toolset on any OS |
| `run_selection.sh` | Full pipeline: rename → tree → FEL/FUBAR/MEME → parse |
| `parse_fel.py` / `parse_fubar.py` / `parse_meme.py` | Convert each HyPhy `.json` into a per-site CSV |

## How to run

```bash
conda env create -f environment.yml
conda activate hyphy
# place pb1_codon_nostop.fasta (in-frame codon alignment, terminal stop removed) here
bash run_selection.sh
```

## Input requirements (why Datamonkey/HyPhy errored before)

1. **In-frame codon alignment** — length a multiple of 3.
2. **Terminal stop codon removed** — 2274 nt → **2271 nt** (757 → ignore stop).
   HyPhy/Datamonkey reject in-frame stop codons. This was the original error.
3. **Clean sequence names** — the pipeline reduces headers like
   `EPI3344803|PB1|A/Colorado/...|...` to just `EPI3344803`. The `/` and `|`
   characters break tree–alignment name matching otherwise.

## Outputs

- `pb1_clean.fasta.FEL.json`, `*.FUBAR.json`, `*.MEME.json` — raw HyPhy results
- `fel_sites.csv`, `fubar_sites.csv`, `meme_sites.csv` — per-site tables
  (site, rates, significance, `is_contact` flag for the 14 interface residues)
- `id_map.tsv` — EPI accession → full strain name

## Key result (this run)

Global tree-wide dN/dS = **0.092** — the gene is overwhelmingly under purifying
selection.

| Method | Purifying (negative) | Positive |
|---|---|---|
| FEL (p ≤ 0.1) | 366 / 757 sites | 2 sites (361, 687) |
| FUBAR (post ≥ 0.9) | 536 / 757 sites | 1 site (687) |
| MEME (p ≤ 0.1, episodic) | — | 3 sites (361, 687, 736) |

**PB1–PB2 contact pocket (14 residues): no positive selection by any method.**
- 7 strongly purifying: 692, 695, 696, 699, 715, 725, 747
- 4 fully invariant (0 substitutions in 500 sequences): 698, 722, 750, 754
  — report descriptively; a site with no substitutions carries no selection signal
- 3 tolerant rim residues: 719, 746, 753
- Positive-selected sites (687, 736) lie in the C-terminal interface region but
  are **not** direct contacts — i.e. the rim adapts, the contact core does not.

## Reproducing on Windows / sharing

- **Windows:** install WSL (`wsl --install`), then follow the same steps inside
  Ubuntu. HyPhy runs cleanly under WSL; a native Windows build is not recommended.
- **Datamonkey (web):** upload `pb1_clean.fasta` (clean EPI names, stop codon
  already removed) and run FEL / FUBAR / MEME there — it builds its own tree.
  The renamed file avoids the header problem seen with raw GISAID names. FUBAR is
  the lightest if the server times out on 500 sequences.
