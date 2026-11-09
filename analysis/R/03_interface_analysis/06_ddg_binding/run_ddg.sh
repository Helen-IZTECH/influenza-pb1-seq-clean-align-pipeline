#!/usr/bin/env bash
# run_ddg.sh -- reproducible setup + binding-ddG scan of one side of the
# PB1-PB2 interface (PDB 3A1G). Works on Linux, WSL and macOS.
#
# Usage:   bash run_ddg.sh PB1      # chain A, 41 fixed residues  -> ddg_all.csv
#          bash run_ddg.sh PB2      # chain B, contact residues   -> ddg_all_PB2.csv
#          bash run_ddg.sh PB2 test # same, but only 2 residues (~15 s) to check setup
#
# Steps:
#   1. Build EvoEF2 from source (once).
#   2. Download + clean 3A1G: keep chains A (PB1) + B (PB2), MSE -> MET.
#   3. RepairStructure (EvoEF2) to fix clashes / missing atoms (once).
#   4. Scan every interface residue x 19 substitutions with compute_ddg.py.
#
# Needs: git, a C++ compiler (macOS: `xcode-select --install`), python3, curl.
set -euo pipefail

PROTEIN="${1:-PB1}"
MODE="${2:-full}"
case "$PROTEIN" in
  PB1) CHAIN=A; OUT=ddg_all.csv ;;
  PB2) CHAIN=B; OUT=ddg_all_PB2.csv ;;
  *) echo "usage: bash run_ddg.sh PB1|PB2 [test]"; exit 1 ;;
esac
EXTRA=""
if [ "$MODE" = "test" ]; then EXTRA="--limit 2"; OUT="test_$OUT"; fi

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
mkdir -p struct

# --- 1. EvoEF2 ------------------------------------------------------------
if [ ! -x EvoEF2/EvoEF2 ]; then
  echo "[1/4] building EvoEF2"
  [ -d EvoEF2 ] || git clone https://github.com/tommyhuangthu/EvoEF2.git
  CXX="${CXX:-g++}"; command -v "$CXX" >/dev/null || CXX=c++
  # gcc wants --fast-math, clang (macOS) wants -ffast-math; -ffast-math works for both
  ( cd EvoEF2 && "$CXX" -O3 -ffast-math -o EvoEF2 src/*.cpp )
fi
EVOEF="$HERE/EvoEF2/EvoEF2"

# --- 2. download + clean structure ---------------------------------------
if [ ! -f struct/3A1G_clean_Repair.pdb ]; then
  [ -f struct/3A1G.pdb ] || curl -fsSL -o struct/3A1G.pdb https://files.rcsb.org/download/3A1G.pdb
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

  # --- 3. repair ----------------------------------------------------------
  echo "[3/4] RepairStructure"
  ( cd struct && "$EVOEF" --command=RepairStructure --pdb=3A1G_clean.pdb >/dev/null )
  [ -f struct/3A1G_clean_Repair.pdb ] || { echo "RepairStructure failed"; exit 1; }
else
  echo "[2-3/4] repaired structure already present, skipping"
fi

# --- 4. ddG scan ----------------------------------------------------------
echo "[4/4] ddG scan of $PROTEIN (chain $CHAIN) -> $OUT"
python3 compute_ddg.py --evoef "$EVOEF" \
    --pdb struct/3A1G_clean_Repair.pdb --chain "$CHAIN" $EXTRA --out "$OUT"

echo "done -> $OUT"
if [ "$PROTEIN" = "PB1" ]; then
  echo "next: Rscript merge_ddg_selection.R   then   Rscript ddg_figure.R"
else
  echo "next: Rscript merge_ddg_selection.R ddg_all_PB2.csv fel_sites_PB2.csv fubar_sites_PB2.csv ddg_vs_selection_PB2.csv"
fi
