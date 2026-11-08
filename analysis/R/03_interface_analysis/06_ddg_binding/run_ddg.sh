#!/usr/bin/env bash
# run_ddg.sh -- one-shot reproducible setup + ddG scan of the PB1-PB2 interface.
#
# Steps:
#   1. Build EvoEF2 from source (once).
#   2. Clean the 3A1G structure: keep chains A (PB1) + B (PB2), MSE -> MET.
#   3. RepairStructure (EvoEF2) to fix clashes / missing atoms.
#   4. Scan all 41 interface residues x 19 substitutions -> ddg_all.csv.
#
# Run from this folder inside WSL/Linux:  bash run_ddg.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
mkdir -p struct

# --- 1. EvoEF2 ------------------------------------------------------------
if [ ! -x EvoEF2/EvoEF2 ]; then
  echo "[1/4] building EvoEF2"
  git clone https://github.com/tommyhuangthu/EvoEF2.git
  ( cd EvoEF2 && g++ -O3 --fast-math -o EvoEF2 src/*.cpp )
fi
EVOEF="$HERE/EvoEF2/EvoEF2"

# --- 2. clean structure ---------------------------------------------------
# Expect raw PDB at struct/3A1G.pdb (download: wget https://files.rcsb.org/download/3A1G.pdb)
echo "[2/4] cleaning structure"
python3 - <<'PY'
keep_ch = {"A","B"}
out = open("struct/3A1G_clean.pdb","w")
for line in open("struct/3A1G.pdb"):
    rec = line[:6].strip()
    if rec in ("ATOM","HETATM"):
        if line[21] not in keep_ch: continue
        resn = line[17:20].strip()
        if rec == "HETATM" and resn != "MSE": continue  # drop waters/ligands
        if resn == "MSE":                               # selenomethionine -> MET
            line = "ATOM  " + line[6:17] + "MET" + line[20:]
            line = line.replace("SE   MET","SD   MET").replace(" SE ", " SD ")
        out.write(line)
    elif rec == "TER":
        out.write(line)
out.write("END\n"); out.close()
print("  wrote struct/3A1G_clean.pdb")
PY

# --- 3. repair ------------------------------------------------------------
echo "[3/4] RepairStructure"
"$EVOEF" --command=RepairStructure --pdb=struct/3A1G_clean.pdb >/dev/null
# EvoEF2 writes 3A1G_clean_Repair.pdb in the cwd
mv -f 3A1G_clean_Repair.pdb struct/ 2>/dev/null || true

# --- 4. ddG scan ----------------------------------------------------------
echo "[4/4] ddG scan (41 x 19 = 779 mutants; takes a while)"
python3 compute_ddg.py --evoef "$EVOEF" \
    --pdb struct/3A1G_clean_Repair.pdb --chain A --out ddg_all.csv

echo "done -> ddg_all.csv"
echo "next: Rscript merge_ddg_selection.R   then   Rscript ddg_figure.R"
