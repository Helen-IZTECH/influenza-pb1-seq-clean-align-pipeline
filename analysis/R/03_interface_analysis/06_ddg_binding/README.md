# 06 · Binding ΔΔG of the PB1–PB2 interface (EvoEF2)

Structural-energetic layer on top of the selection analysis (stage 05).
For every PB1 interface contact residue we compute how much each possible
mutation destabilises PB1–PB2 binding, then ask whether evolutionary
conservation tracks binding energetics — and show, with Sugiyama et al.
(2009) experimental data, why it only partly does.

## Biological question

A conserved residue is not automatically a resistance-robust drug-target
residue (PA-I38 is conserved yet baloxavir resistance I38T emerges). The
sharper question: **are the escape mutations functionally non-viable?**
Binding ΔΔG tests one dimension of that — does mutating the residue break
the complex?

## Pipeline

Works on Linux, WSL and **macOS** (needs `git`, a C++ compiler, `python3`,
`curl`, `R`; on a Mac run `xcode-select --install` once).

```bash
# 0. quick check that your setup works (2 residues, ~15 s)
bash run_ddg.sh PB2 test

# 1. PB1 side: build EvoEF2, clean + repair 3A1G, scan 41 x 19 = 779 mutants
bash run_ddg.sh PB1                 # -> ddg_all.csv

# 1b. PB2 side: contact residues are found automatically (4.5 A to PB1)
bash run_ddg.sh PB2                 # -> ddg_all_PB2.csv   (~3 min, 25 x 19 = 475 mutants)

# 2. collapse to per-residue, join FEL/FUBAR selection, classify
Rscript merge_ddg_selection.R       # PB1 -> ddg_vs_selection.csv
Rscript merge_ddg_selection.R ddg_all_PB2.csv fel_sites_PB2.csv \
        fubar_sites_PB2.csv ddg_vs_selection_PB2.csv        # PB2

# 3. figures (PB1)
Rscript ddg_figure.R                # -> FIG_ddg.png (ddG vs conservation)
Rscript convergence_figure.R        # -> FIG_convergence.png (3 lines of evidence)
```

Inputs expected in this folder: the selection tables `fel_sites.csv`,
`fubar_sites.csv` from stage 05 (PB2: `fel_sites_PB2.csv`,
`fubar_sites_PB2.csv`, with PB2 residue numbers in the `site` column).
The PDB is downloaded automatically; if your network blocks it, save
https://files.rcsb.org/download/3A1G.pdb as `struct/3A1G.pdb` yourself.

### Troubleshooting (macOS)

- `g++: command not found` -> run `xcode-select --install`.
- EvoEF2 compile warnings about `printf` format are harmless; only an
  `error:` line matters (send it to the maintainer).
- Reproducibility check: `bash run_ddg.sh PB1 test` should give
  `685,D,A,-0.52` as the first data row of `test_ddg_all.csv`.
- PB2 residue numbers <= 0 in 3A1G (expression tag) are skipped automatically.

## Method notes

- **Structure**: PDB 3A1G (A/Puerto Rico/8/34), chain A = PB1 (685–757),
  chain B = PB2. Selenomethionine (MSE) → MET; `RepairStructure` before scanning.
- **ΔΔG_bind** = E_bind(mutant) − E_bind(wild type). Positive = destabilising.
  WT binding dG ≈ −162.8.
- **Mutant syntax** (EvoEF2): `{WT}{Chain}{Pos}{Mut}` e.g. `VA715S`.
- **Per residue** we keep `mean_ddG` (average over 19 substitutions) and
  `min_ddG` (the least-destabilising substitution = the *best escape* at the
  binding level).
- **Classes** (from stage 05): `invariant` (FEL α=β=0), `purifying`
  (FUBAR posterior ≥ 0.9), `tolerant` (the rest).

## Results (honest summary)

Per-class mean ΔΔG: tolerant **1.9** < invariant **3.0** < purifying **3.4**.
Spearman(mean_ddG, FUBAR purifying) = **0.19** (weak); Spearman(mean_ddG, β)
= **−0.37** (moderate — residues evolution lets vary are less destabilising).

The correlation is weak **by design, not by failure**: interface residues are
constrained on several independent axes (binding, fold, catalysis, inter-subunit
dynamics), and binding ΔΔG sees only one. The four residues Sugiyama tested as
single mutants make this concrete (`FIG_convergence.png`):

| residue | nature | ΔΔG min | Sugiyama binding | Sugiyama function | lock is… |
|---------|--------|---------|------------------|-------------------|----------|
| L695    | conserved | +1.3 | lost     | OK       | binding/structural |
| F699    | conserved | +1.6 | weak     | higher   | binding/structural |
| I750    | conserved | +2.7 | lost     | OK       | binding/structural |
| **V715**| conserved | **−2.3** | **retained** | **−80% (dead)** | **dynamics/signalling** |

**V715 is the key case**: a binding-viable escape exists (ΔΔG min < 0;
Sugiyama confirms binding is retained), yet evolution still locks the residue —
because that escape destroys polymerase function through a mechanism binding
energetics cannot see. This is the essence of resistance-robustness: the
interface is multiply constrained, so viable escape must satisfy *all* axes at
once, and evolution (cross-validated by experiment) reads that integral better
than any single computational metric.

## Files

| file | what |
|------|------|
| `run_ddg.sh` | `PB1` or `PB2` [`test`]: build EvoEF2, fetch/clean/repair 3A1G, run the scan |
| `compute_ddg.py` | per-mutant EvoEF2 BuildMutant + ComputeBinding driver (`--chain A|B`, `--residues default|auto|list`) |
| `merge_ddg_selection.R` | per-residue ΔΔG + FEL/FUBAR join + classification |
| `ddg_figure.R` | ΔΔG vs conservation (box by class + scatter) |
| `convergence_figure.R` | evolution × structure × Sugiyama experiment matrix |
| `../../../results/ddg_all.csv` | 779 rows: pos, wt, mut, ddG_bind |
| `../../../results/ddg_vs_selection.csv` | per-residue: pos, class, beta, purifying_post, mean_ddG, min_ddG |

## Caveats

- PB1 uses the fixed 41-residue list from the paper. `--residues auto` (4.5 A
  heavy-atom contact) gives 39 PB1 residues, not the identical set (it misses
  741, 745, 749, 752 and adds 729, 735), so the fixed list is kept for PB1.
  PB2 uses `auto`.
- EvoEF2 ΔΔG is a **prediction**, not a measurement; Sugiyama's assays are the
  experimental ground truth where the two disagree.
- `mean_ddG` averages over all 19 substitutions and over-weights unlikely bulky
  changes; for the *escape* question `min_ddG` is the more meaningful statistic.
- Single-residue binding ΔΔG does not capture fold stability, catalysis, or the
  in-vivo fitness that evolutionary conservation integrates.
