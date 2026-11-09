#!/usr/bin/env bash
# Site-wise dN/dS selection analysis of a codon alignment (HyPhy). PB1 (default) or PB2.
# Reproduces FEL, FUBAR and MEME on the 500-sequence representative set.
# Run inside the conda env defined in environment.yml:
#   conda env create -f environment.yml
#   conda activate hyphy
#   bash run_selection.sh          # PB1  (needs pb1_codon_nostop.fasta)
#   bash run_selection.sh PB2      # PB2  (needs pb2_codon_nostop.fasta)
set -euo pipefail

PROTEIN="${1:-PB1}"
case "$PROTEIN" in PB1|PB2) ;; *) echo "usage: bash run_selection.sh [PB1|PB2]"; exit 1 ;; esac
LOW=$(echo "$PROTEIN" | tr 'A-Z' 'a-z')

RAW=${LOW}_codon_nostop.fasta  # in-frame codon alignment, terminal stop codon removed
ALN=${LOW}_clean.fasta
TREE=${LOW}_clean.treefile
MAP=id_map.tsv; [ "$PROTEIN" = "PB2" ] && MAP=id_map_pb2.tsv

# 1. Simplify sequence names to the EPI accession (first |-field) so the
#    alignment names and the tree tip labels match (HyPhy requires this).
awk '/^>/{split(substr($0,2),a,"|"); print ">"a[1]; next} {print}' "$RAW" > "$ALN"
grep '^>' "$RAW" | sed 's/^>//' | awk -F'|' '{print $1"\t"$0}' > "$MAP"

# 2. Maximum-likelihood tree (GTR+G).
iqtree -s "$ALN" -m GTR+G -T AUTO --prefix ${LOW}_clean -redo

# 3. Selection analyses.
hyphy fel   --alignment "$ALN" --tree "$TREE"   # pervasive selection, frequentist
hyphy fubar --alignment "$ALN" --tree "$TREE"   # pervasive selection, Bayesian
hyphy meme  --alignment "$ALN" --tree "$TREE"   # episodic positive selection

# 4. Parse JSON outputs into per-site CSV tables.
python3 parse_fel.py "$PROTEIN"
python3 parse_fubar.py "$PROTEIN"
python3 parse_meme.py "$PROTEIN"
