#!/usr/bin/env bash
# Site-wise dN/dS selection analysis of the PB1 codon alignment (HyPhy).
# Reproduces FEL, FUBAR and MEME on the 500-sequence representative set.
# Run inside the conda env defined in environment.yml:
#   conda env create -f environment.yml
#   conda activate hyphy
#   bash run_selection.sh
set -euo pipefail

RAW=pb1_codon_nostop.fasta     # in-frame codon alignment, terminal stop codon removed
ALN=pb1_clean.fasta
TREE=pb1_clean.treefile

# 1. Simplify sequence names to the EPI accession (first |-field) so the
#    alignment names and the tree tip labels match (HyPhy requires this).
awk '/^>/{split(substr($0,2),a,"|"); print ">"a[1]; next} {print}' "$RAW" > "$ALN"
grep '^>' "$RAW" | sed 's/^>//' | awk -F'|' '{print $1"\t"$0}' > id_map.tsv

# 2. Maximum-likelihood tree (GTR+G).
iqtree -s "$ALN" -m GTR+G -T AUTO --prefix pb1_clean -redo

# 3. Selection analyses.
hyphy fel   --alignment "$ALN" --tree "$TREE"   # pervasive selection, frequentist
hyphy fubar --alignment "$ALN" --tree "$TREE"   # pervasive selection, Bayesian
hyphy meme  --alignment "$ALN" --tree "$TREE"   # episodic positive selection

# 4. Parse JSON outputs into per-site CSV tables.
python3 parse_fel.py
python3 parse_fubar.py
python3 parse_meme.py
